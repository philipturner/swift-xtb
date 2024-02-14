#!/bin/bash

# Whether the next argument is the user.
# Example: "philipturner" in "/Users/philipturner/miniforge3/lib"
next_argument_user=false

# Whether the user has been specified.
user_specified=false

# The name of the user.
user_name=""

# Whether the next argument is the filter.
next_argument_filter=false

# Whether the filter has been specified.
filter_specified=false

# The flags for the filter.
filter_args=""

# Check whether the arguments parse correctly.
if [[ $# -gt 4 ]]; then
  echo "Too many arguments."
  invalid_input=true
elif [[ $# != 0 ]]; then
  invalid_input=false
 
  if [[ $invalid_input == false ]]; then
    for param in "$@"; do
      if [[ $param == "--user" ]]; then
        if [[ $next_argument_user == true ]]; then
          echo "Duplicate argument '--user'."
          invalid_input=true
        elif [[ $user_specified == true ]]; then
          echo "Duplicate argument '--user'."
          invalid_input=true
        elif [[ $next_argument_filter == true ]]; then
          echo "Invalid use of '--filter'."
          invalid_input=true
        else
          next_argument_user=true
        fi
      elif [[ $param == "--filter" ]]; then
        if [[ $next_argument_filter == true ]]; then
          echo "Duplicate argument '--filter'."
          invalid_input=true
        elif [[ $filter_specified == true ]]; then
          echo "Duplicate argument '--filter'."
          invalid_input=true
        elif [[ $next_argument_user == true ]]; then
          echo "Invalid use of '--user'."
          invalid_input=true
        else
          next_argument_filter=true
        fi
      elif [[ $next_argument_user == true ]]; then
        next_argument_user=false
        user_specified=true
        user_name="${param}"
      elif [[ $next_argument_filter == true ]]; then
        next_argument_filter=false
        filter_specified=true
        filter_args="${param}"
      else
        echo "Unrecognized argument '${param}'."
        invalid_input=true
      fi
    done
  fi
else
  echo "No arguments found."
  invalid_input=true
fi

# If the arguments parse correctly, check whether the path is correct.
if [[ $invalid_input == false ]]; then
  if [[ $user_specified == true ]]; then
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
      # Don't know where the Linux path for libxc is yet.
      export DFTD4_LIBRARY_PATH="/home/${user_name}/miniconda3/lib"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
      export XC_LIBRARY_PATH="/opt/homebrew/lib"
      export DFTD4_LIBRARY_PATH="/Users/${user_name}/miniforge3/lib"
    else
      echo "Unrecognized OS for argument '--user'."
      exit -1
    fi
  else
    echo "No user specified."
    invalid_input=true
  fi
fi

# Return early if the arguments are incorrect.
if [[ $invalid_input == true ]]; then
  echo "Usage: test.sh [--user USER] [--filter FILTER]"
  exit -1
fi

# Run the Swift package tests.
if [[ $filter_specified == true ]]; then
  swift test -Xswiftc -Ounchecked --filter "$filter_args"
else
  swift test -Xswiftc -Ounchecked
fi
