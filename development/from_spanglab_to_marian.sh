#!/bin/bash

PATH_TO="/Users/marian/Documents/"
SOURCE="${PATH_TO}/git/DTD_spanglab"
TARGET="${PATH_TO}/git/DTD"

rsync -av ${SOURCE}/ ${TARGET}/ --exclude-from=exclude_from_rsync
