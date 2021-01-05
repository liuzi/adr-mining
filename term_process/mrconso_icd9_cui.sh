#!/usr/bin/env bash
folder=$1

pushd $folder
grep "|ENG|" MRCONSO.RRF | grep "ICD9CM" > MRCONSO_ICD9.RRF
popd