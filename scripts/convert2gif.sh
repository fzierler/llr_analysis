name=$1

convert -alpha deactivate -delay 9 -loop 0 -density 100 $name.pdf $name-tmp.gif
gifsicle --colors 256 -i $name-tmp.gif -O3 -o $name.gif
rm $name-tmp.gif