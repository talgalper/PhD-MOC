#!/bin/bash

# fix ubuntu/linux detection and install

### Checking and installing fpocket requirements ###

# Detect the user's system (Mac or Linux)
if [[ "$(uname)" == "Darwin" ]]; then
    echo "Mac detected"
    echo "Setting up fpocket for Mac..."

    # Check if MacPorts is installed
    if [ -x "$(command -v port)" ]; then
      echo "MacPorts is installed"
    else
      echo "Must install MacPorts from: https://www.macports.org/"
    fi

    # Check if the libnetcdf-dev package is installed
    dpkg -s libnetcdf-dev &> /dev/null
    if [ $? -eq 0 ]; then
        echo "The libnetcdf-dev package is already installed."
    else
        sudo port install netcdf
        export LIBRARY_PATH=/opt/local/lib
        if [ $? -eq 0 ]; then
            echo "The libnetcdf-dev package was successfully installed."
        else
            echo "There was an error installing the libnetcdf-dev package."
        fi
    fi

    # Check if the fpocket repository has already been cloned
    if [ -d "fpocket" ]; then
      echo "The fpocket repository has already been cloned."
    else
      # Clone the fpocket repository
      git clone https://github.com/Discngine/fpocket.git
      if [ $? -eq 0 ]; then
        echo "The fpocket repository was successfully cloned."
      else
        echo "There was an error cloning the fpocket repository."
      fi
    fi

    # Compile fpocket
    cd "fpocket"
    if [ -f "bin/fpocket" ]; then
      echo "fpocket has already been compiled."
    else
      make ARCH=MACOSXX86_64
      if [ $? -eq 0 ]; then
        echo "fpocket was successfully compiled."
      else
        echo "There was an error compiling fpocket."
        exit 1
      fi
      sudo make install
      if [ $? -eq 0 ]; then
        echo "fpocket was successfully installed."
      else
        echo "There was an error installing fpocket."
        exit 1
      fi
    fi

  # Test fpocket
  pytest
  if [ $? -eq 0 ]; then
    echo "fpocket passed all tests."
  else
    echo "fpocket failed some tests."
  fi


elif [[ "$(lsb_release -si)" == "Ubuntu" ]]; then
    echo "Ubuntu detected"
    echo "Setting up fpocket for Ubuntu..."

    dpkg -s libnetcdf-dev &> /dev/null
    if [ $? -eq 0 ]; then
        echo "The libnetcdf-dev package is already installed."
    else
        # Install the libnetcdf-dev package
        sudo apt-get install -y libnetcdf-dev
        if [ $? -eq 0 ]; then
            echo "The libnetcdf-dev package was successfully installed."
        else
            echo "There was an error installing the libnetcdf-dev package."
        fi
    fi

    # Check if the fpocket repository has already been cloned
    if [ -d "fpocket" ]; then
      echo "The fpocket repository has already been cloned."
    else
      # Clone the fpocket repository
      git clone https://github.com/Discngine/fpocket.git
      if [ $? -eq 0 ]; then
        echo "The fpocket repository was successfully cloned."
      else
        echo "There was an error cloning the fpocket repository."
      fi
    fi

    # Compile fpocket
    cd "fpocket"
    if [ -f "bin/fpocket" ]; then
      echo "fpocket has already been compiled."
    else
      make
      if [ $? -eq 0 ]; then
        echo "fpocket was successfully compiled."
      else
        echo "There was an error compiling fpocket."
        exit 1
      fi
      sudo make install
      if [ $? -eq 0 ]; then
        echo "fpocket was successfully installed."
      else
        echo "There was an error installing fpocket."
        exit 1
      fi
    fi

else
  echo "Unknown system detected. Exiting..."
  exit 1
fi

echo "fpocket setup complete"
