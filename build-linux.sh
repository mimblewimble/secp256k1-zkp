#!/bin/bash

rm -rf "build"
cmake -H. -Bbuild
cmake --build ./build --target secp256k1 --config Release