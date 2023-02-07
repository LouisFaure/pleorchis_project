gget blast -l 1 -csv $(awk '/^>/{if(NR==1){print}else{printf("\n%s\n",$0)}next} {printf("%s",$0)} END{printf("\n")}' $1 | sed -n '2p' | less) > $2
