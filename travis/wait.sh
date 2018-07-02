#!/bin/bash
if [[ $test_complete == false ]]; then
     sleep 2700;
     echo "slept for 45 min";
     echo "test is complete?" $test_complete;
     echo "tail output.txt";
     tail -n 100 output.txt;
fi