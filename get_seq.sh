awk -v seq=$1 -v RS='>' '$1 == seq {{print RS $0}}' $2
