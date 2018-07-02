#!/bin/bash
if [ ! -f done.txt ]; then
  for ((i=0; i < 45; i++)); do
    sleep 60
    printf "…"${i}
  done
  echo "slept for 45 min";
  echo "tail output.txt";
  tail -n 100 output.txt;
fi