#!/bin/bash
if [[ $test_complete == false ]]; then
  for ((i=0; i < 45; i++)); do
    sleep 60
    printf "…"${i}
  done
  echo "slept for 45 min";
  echo "test is complete?" $test_complete;
  echo "tail output.txt";
  tail -n 100 output.txt;
fi