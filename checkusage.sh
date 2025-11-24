bash

#Bash script to monitor changes to file
# THIS SCRIPT HAS NO DEFAULT PARAMETERS AND SHOULD NOT BE RAN INDEPENDENT OF SpecifiedBashScriptForTimeout.sh
# We create a function that will check every 5 minutes if a filesize has been changed or not
# If the file has not been changed we kill the pid related to the programme creating that file
# As we assume the process has either A finished or B is hanging 
#!/bin/bash
checkUsage() # Note that U is in caps
{
while true
do
        sleep 300
        fileSize=$(stat -c%s $1)
        sleep 240;
        fileSizeNew=$(stat -c%s $1)

        if [ "$fileSize" == "$fileSizeNew" ]
        then
           echo -e  "[Notice : ] no changes noted in $1 : gracefully exiting"
           kill -HUP $2
           return 1 # leave the while loop
        fi

done
}
checkUsage $1 $2   # this function takes in arguements of the form $1 /some/file/we/are/watching $2 pid of the function creating $1