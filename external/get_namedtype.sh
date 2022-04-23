#!/bin/sh

git clone --branch master --depth 1 https://github.com/joboccara/NamedType.git NamedType    
cd NamedType
rm -rf test cmake .git .gitignore .clang-format .travis.yml 
