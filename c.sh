#!/bin/bash
javac -cp *:. rt/Main.java
if [ $? -eq 0 ]
then
mkdir output
java -ea -cp *:. rt.Main
fi

