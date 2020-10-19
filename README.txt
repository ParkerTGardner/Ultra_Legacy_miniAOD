To compile (check for compilation errors):

make clean; make


To run a single job: ./bin/MyClass.exe <List of Input Files> <Job Number> <Total Jobs>
Example for testing: ./bin/MyClass.exe fileLists/UL_D_Files_ak8_test.txt 0 1

To run parallel: ./bash/runMyClass.sh <List of Input Files> <Total jobs>
Output and logs will appear in unmergedOutputs directory
Try to keep Total Jobs less than the number of Cores on a machine (16?).


To kill running parallel jobs (which run in background)
top
kill <process ID number>
y
