#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
echo $0
echo $1
javac ${DIR}TangleHtml.java
javac ${DIR}Tangle.java

if [ $? -eq 0 ] 
then
java -cp ${DIR} TangleHtml $1
java -cp ${DIR} Tangle $1
fi

