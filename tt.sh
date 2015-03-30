#!/bin/bash
./t.sh all.java
if [ $? -eq 0 ] 
then
./c.sh
fi

