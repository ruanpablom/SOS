#!/bin/bash

[ -z "$1" ] && echo "Forneça o nome de um arquivo como parâmetro" && echo && \
echo "./opComValorEmArq.sh nomedoarquivo" && echo && exit

! [ -f "$1" ] && echo "O parâmetro [ $1 ] não é um arquivo" && exit

SOMA=0

# Conta as linhas que tenham conteúdo
LINHAS=$(cat "$1" | egrep '[^(^$)]' | egrep '[0-9]+' | wc -l | cut -d" " -f 1)

for i in $(seq $LINHAS)
do
   VALORDALINHA=$(cat "$1" | egrep -o '[0-9]+' | head -n $i | tail -n 1)
   #soma
   ! [ -z "$VALORDALINHA" ] && SOMA=$(($SOMA+$VALORDALINHA))
done

#media
MEDIA=$(($SOMA/$LINHAS))
echo "$SOMA"
echo $MEDIA >> mediasCORES.txt;
#echo  >> mediasCORES.txt;