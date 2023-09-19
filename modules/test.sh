#!/bin/bash


a=10
b=20


test=$(awk -v var1=3 -v var2=4 'BEGIN { print  ( var1 / var2 ) }')


echo $test
