for name in $@
do
    extension=${name#*.}
    short_name=${name%%_*}
    mv $name $short_name-2.$extension
done
