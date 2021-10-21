#!/bin/bash
if [ -z ${ELECTRON_SPECIFIC_HEAT_DIR+x} ] ; then
	SOURCE="${BASH_SOURCE[0]}"
	while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
	  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
	  SOURCE="$(readlink "$SOURCE")"
	  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
	done
	export ELECTRON_SPECIFIC_HEAT_DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
	export PATH+=:$ELECTRON_SPECIFIC_HEAT_DIR/bin:$ELECTRON_SPECIFIC_HEAT_DIR/script:$ELECTRON_SPECIFIC_HEAT_DIR/tensorflow
	if [ $(uname) = "Darwin" ]; then
		export DYLD_LIBRARY_PATH+=:$ELECTRON_SPECIFIC_HEAT_DIR/lib
	else
		export LD_LIBRARY_PATH+=:$ELECTRON_SPECIFIC_HEAT_DIR/bin
	fi
	alias eshcd='cd $ELECTRON_SPECIFIC_HEAT_DIR'
	alias eshmake='cd $ELECTRON_SPECIFIC_HEAT_DIR/build; make -j4; cd -'
fi