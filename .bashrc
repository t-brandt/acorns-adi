notify=
MAILCHECK=120
history_control=ignoredups
no_exit_on_failed_exec=
ignoreeof=1
export PRINTER=ps
#
# Check for an interactive shell
#

#. /usr/peyton/common/licensed/idl64/bin/idl_setup.bash

# Set python path to point to scipy, numpy, pyfits, etc.
# Wild cards make this *slightly* flexible to a different path for 
#   a new version
for D in /usr/peyton/utils/Python*/lib/python*/site-packages ; do
    export PYTHONPATH=$D
done
export PYTHONPATH=$PYTHONPATH:~/Dropbox/SEEDS/destripe_tim/destripe

# Search /usr/peyton/utils for PATH directories first
# Also grab library and manual paths.
# I have not been able to make include work well here; I just added
#   a flag as an alias to my gcc

export PATH=""
for D in /usr/peyton/utils/* ; do
    if [ -d $D/bin ]; then
       [[ "$PATH" =~ "(^|:)$D/bin($|:)" ]] || \
       export PATH=${PATH:+$PATH:}$D/bin
    fi

    if [ -d $D/share/man ]; then
       [[ "$MANPATH" =~ "(^|:)$D/share/man($|:)" ]] || \
       export MANPATH=${MANPATH:+$MANPATH:}$D/share/man
    fi

    if [ -d $D/man ]; then
       [[ "$MANPATH" =~ "(^|:)$D/man($|:)" ]] || \
       export MANPATH=${MANPATH:+$MANPATH}$D/man:
    fi

    if [ -d $D/lib ]; then
       [[ "$LPATH" =~ "(^|:)$D/lib($|:)" ]] || \
       export LPATH=${LPATH:+$LPATH}$D/lib:
    fi

    if [ -d $D/lib ]; then
       [[ "$LD_LIBRARY_PATH" =~ "(^|:)$D/lib($|:)" ]] || \
       export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH}$D/lib:
    fi
done

export LPATH=${LPATH}:/usr/peyton/utils/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/peyton/utils/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/peyton/utils/ATLAS/Linux_UNKNOWNSSE2_4/lib

# Add these to PATH last, so that we search for a special installation 
#   in /usr/peyton/utils first

export PATH=$PATH:/usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin:/usr/peyton/bin:/usr/X11R6/bin:/opt/SUNWspro/bin:/usr/lang:/usr/ucb:/scr/wakusei1/software/ql2:/u/tbrandt/bin

#export JAVA_HOME=/usr/peyton/jre1.6.0_05
#export JAVA_HOME=/u/tbrandt/j/jdk1.6.0_14
export DUST_DIR=/u/schlegel/dustpub
export set CVS_RSH=ssh
export set CVSROOT=:ext:brandt@battra.lbl.gov:/home/ccse_src/cvsroot
export set CVSEDITOR=emacs

# Automatically link math library and include the utils files with gcc
#alias gcc='gcc -I/usr/peyton/utils/include -I/usr/peyton/utils/ATLAS/include -lm'
export C_INCLUDE_PATH=${C_INCLUDE_PATH}:/usr/peyton/utils/include
export C_INCLUDE_PATH=${C_INCLUDE_PATH}:/usr/peyton/utils/ATLAS/include

# Colors for ls command
export LS_OPTIONS='--color=auto'
eval `dircolors`
alias ls='ls $LS_OPTIONS'

if [ "$PS1" ]; then
  CDPATH=.:~:~/sm
  HISTSIZE=80
  set_prompt () {
	PS1=`hostname|colrm 3 99`:`basename $PWD`">"
  }
  set_prompt

  cd () {
	old=$PWD; builtin cd $1;
	set_prompt
  }
  back () {
	local backd=$old; old=$PWD; builtin cd $backd;
	set_prompt
  }
  pd () {
	if [ "$1" = "" ]; then
		pushd
	else
		_old=$old; pushd $1;
	fi
	set_prompt
  }
  popd () {
	old=$_old; unset _old
	builtin popd
	set_prompt
  }
fi
if [ -f ~/.bash_aliases ]; then
	source ~/.bash_aliases
fi

. ~condor/condor-setup.sh

# start Dropbox
if ps ax | grep -v grep | grep dropbox > /dev/null
then
    echo "" > /dev/null
else
    echo "Starting Dropbox..."
    ~/.dropbox-dist/dropboxd &
fi

export IDL_PATH=${IDL_PATH}:+/scr/wakusei1/users/tbrandt/adi
export IDL_PATH=${IDL_PATH}:/usr/linux/common/peyton/licensed/idl71/lib 
export IDL_PATH=${IDL_PATH}:/u/schlegel/dustpub/CodeIDL #Schlegel's dust maps
export IDL_PATH=${IDL_PATH}:/u/tbrandt/idl #my stuff
export IDLSPEC2D_DIR=/u/schlegel/idlspec2d

#setup sloan stuff
source /u/dss/products/eups/bin/setups.sh
setup -v idlutils
#setup -v photoop v1_9_11
setup -v idl v7_1

