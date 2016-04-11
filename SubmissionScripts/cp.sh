function joinStrings { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }
#cp " joinStrings '\ ' $@ " /home/heath/FastQC2/ 
cp "$@" /home/heath/FastQC2/ 

