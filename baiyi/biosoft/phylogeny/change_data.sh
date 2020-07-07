#!/bin/bash -
#===============================================================================
#
#          FILE: change_data.sh
#
#         USAGE: ./change_data.sh
#
#   DESCRIPTION: 
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Yuan Hanyu (Bioinformatics Engineer), yhy.119@foxmail.com
#  ORGANIZATION: Beijing Computing Center
#       CREATED: 05/27/2019 02:56:52 PM
#      REVISION:  ---
#===============================================================================

set -o nounset                                  # Treat unset variables as an error

ls | sed -n '/^C/p' |  awk '{print "mv "$0" N"$0}' | bash

