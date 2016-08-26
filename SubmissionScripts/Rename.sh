for name in $@
do
    extension=${name#*.}
    short_name=${name%%_*}
    echo "mv $name $short_name.$extension"
done
