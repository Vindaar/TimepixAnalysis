#!/usr/bin/env bash

# NOTE: this script *will* possibly modify:
# - $HOME/src/nim_git_repo for the Nim compiler (name such that it won't conflict
#   with other dirs by chance
# - append a line in .bashrc / .zsh to export PATH incl. Nim compiler

# The script runs using `bash`, but adds the Nim compiler to the PATH of the
# shell defined by the `SHELL` variable

# check if Nim in PATH
if hash nim 2>/dev/null; then
    echo "Found Nim compiler"
else
    mkdir -p $HOME/src
    cd $HOME/src
    git clone https://github.com/nim-lang/nim nim_git_repo
    cd nim_git_repo
    sh build_all.sh
    bin/nim c koch
    ./koch boot -d:release
    ./koch tools

    toApp='export PATH=$PATH:$HOME/src/nim_git_repo/bin'
    if [[ $SHELL == *"bash"* ]]; then
        echo $toApp | tee $HOME/.bashrc
        source $HOME/.bashrc
    elif [[ $SHELL == *"zsh"* ]]; then
        echo $toApp | tee $HOME/.zshrc
        source $HOME/.zshrc
    else
        echo "Neither bash nor zsh used as SHELL. Append "
        echo "$HOME/src/nim_git_repo/bin"
        echo "to your PATH manually!"
    fi
fi

# sanity check for compiler again to make sure above worked
if hash nim 2>/dev/null; then
    echo "Building TimepixAnalysis!"
    nim e buildTpa.nims
fi
