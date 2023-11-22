#!/bin/sh

git clone --branch master --depth 1 https://github.com/carmeloevoli/sophia_next.git sophianext    
cd sophianext
rm -rf test cmake .git .gitignore .clang-format .travis.yml 
