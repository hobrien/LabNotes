function joinStrings { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }
#cp " joinStrings '\ ' $@ " /home/mpmho/FastQC2/ 
cp "$@" /home/mpmho/FastQC2/ 

