PROC_DIR=$1
TMP_DIR=$2
for f in "$PROC_DIR"/out_*;
do
  filename="${f##*/}"
  newfilename=${filename//"out"/}
  grep "\[WILSONFLOW\]" "$PROC_DIR"/$filename | cut -d"=" -f2 | awk '{if($1==0.0) i++; print i,$0;}' > "$TMP_DIR"/WF"$newfilename"	
done
