#! /usr/bin/env bash

# this removes the :*_latest tag and replaces with versioned container
  for i in harpy/snakefiles/*.smk; do
    sed -i "s/_latest/_${1}/g" $i
  done

# this removes the :latest tag and replaces with versioned container
sed -i "s/-dev//g" harpy/common/version.py
sed -i "s/-dev//g" pyproject.toml     
