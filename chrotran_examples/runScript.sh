#!/bin/bash

date
pwd
#

for i in "$@"
do
case $i in
    -t=*|--testdir=*)
    TESTDIR="${i#*=}"
    ;;
    -f=*|--filename=*)
    FILENAME="${i#*=}"
    ;;
    -np=*|--nproc=*)
    NPROC="${i#*=}"
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
done


MYDIR=/lclscratch/sach/Programs/pflotran-chrome-bio/chrome-bio-tests/
INPUT_FILE=${FILENAME}.in
PFLOTRAN_EXE=/lclscratch/sach/Programs/pflotran-chrome-bio/src/pflotran/pflotran
OUTPUT_FILE=${FILENAME}.log
NPROCS=16

cd ${MYDIR}/${TESTDIR}

if [ -f "${OUTPUT_FILE}" ]; then
  echo "log file exists!"
  echo "deleting that file"
  rm -rf ${OUTPUT_FILE}
fi

if [ "${NPROC}" -gt 1 ]
then
echo "mpirun -np ${NPROCS} ${PFLOTRAN_EXE} -pflotranin ${INPUT_FILE} >& ${OUTPUT_FILE}"
mpirun -np ${NPROCS} ${PFLOTRAN_EXE} -pflotranin ${INPUT_FILE} >& ${OUTPUT_FILE}
else
echo "${PFLOTRAN_EXE} -pflotranin ${INPUT_FILE} >& ${OUTPUT_FILE}"
${PFLOTRAN_EXE} -pflotranin ${INPUT_FILE} >& ${OUTPUT_FILE}
fi


