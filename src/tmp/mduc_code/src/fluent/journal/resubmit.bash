#!/bin/bash

# ----- arguments
case_file=$1
current_n=$2
current_t=`printf "%.0f" "$3"`
event=$4

echo " ====================================================================== "
echo " "
echo " Resubmitting case from: $case_file"
echo " Current time step:      $current_n"
echo " Current time:           $current_t"
if [ -n "$event" ]; then
    echo " From dynamic event:     $event"
fi
echo " "

# ----- get the case parameters from the scheme file
cwd=`pwd`
case_name=`grep 'case-name' $case_file | awk '{print substr($3,2)}'`
case_name=${case_name%\")*}
case_map=`grep 'case-map' $case_file | awk '{print substr($3,2)}'`
case_map=${case_map%\")*}
case_spec=`grep '(define case-spec' $case_file | awk '{print substr($3,2)}'`
case_spec=${case_spec%\")*}
case_id=`grep 'case-id' $case_file | awk '{print substr($3,2)}'`
case_id=${case_id%\")*}
case_dir="$case_name/$case_map/$case_spec"
mduc_file=`grep 'mduc-input' $case_file | awk '{print substr($3,2)}'`
mduc_file=${mduc_file%\")*}

# get the event file and change it if necessary
event_file=`grep 'event-file' $case_file | awk '{print substr($3,2)}'`
event_file=${event_file%\")*}
if [ -n "$event" ]; then
    # turn off the resubmit action
    new_events=`printf "%s_rs%05i.scm" "${event_file%.scm}" "$current_n"`
    sed -e '/'$event'/c;; !!!!! event removed: '$event' !!!!!' \
        < $case_spec/init/$event_file > $case_spec/init/$new_events

    # reassign the event file
    event_file=$new_events
fi

echo " Case info:"
echo "   name:            $case_name"
echo "   test map:        $case_map"
echo "   operating point: $case_spec"
echo "   run id:          $case_id"
echo "   events:          $event_file"
echo " "

# reconstruct the file prefix (be sure to change this if necessary
file_prefix="$case_name-$case_map_$case_spec-$case_id"

# ----- create a new mduc file
# define the new MDUC restart file
mduc_restart=`printf "%s_%05i" "$file_prefix" "$current_n"`
mduc_restart="$cwd/$case_spec/rd_$case_id/mduc/$mduc_restart"
mduc_file_n=`printf "mduc_rs%05i.inp" "$current_n"`
sed -e '/read-data-file/cread-data-file '$mduc_restart'.rif' \
    -e '/read-flow-file/cread-flow-file '$mduc_restart'.flo' \
    < $case_spec/init/$mduc_file > $case_spec/init/$mduc_file_n

echo " MDUC restart file copied: $mduc_file --> $mduc_file_n"
echo " "

# ----- create a new input file with the correct restart time
# define the name of the scheme file
scm_restart=`printf "%s_rs%05i.scm" "${case_file%.scm}" "$current_n" `

# change the restart parameters and the start time
sed -e '/case-restart/c(define case-restart #t)' \
    -e '/mduc-restart/c(define mduc-restart #t)' \
    -e '/mduc-input/c(define mduc-input \"'$mduc_file_n'\")' \
    -e '/event-file/c(define event-file \"'$event_file'\")' \
    -e '/nt-start/c\(define nt-start '${current_t//[A-Za-z]/}')' \
    < $case_file > $scm_restart

echo " Scheme file created: $case_file --> $scm_restart"
echo " "

# ----- create a new queue submission script
# create the new file name
qsub_file_n=`printf "sub_%s_rs%05i.sge" "$case_spec" "$current_n"`

# modify the submission script from the template
sed -e '/INPUT_FILE=/cINPUT_FILE=\"'$scm_restart'\"' \
    -e 's/XXX/'${case_file%.scm}_rs$current_n'/' \
    < sub.sge > $qsub_file_n

echo " Queue submission script created: sub.sge --> $qsub_file_n"
echo " "

# ----- resubmit the job
chmod +x $qsub_file_n
qsub $qsub_file_n

echo " "
echo " ====================================================================== "
# ===== end file ===== #
