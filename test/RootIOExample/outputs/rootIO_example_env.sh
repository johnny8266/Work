if [ -z "${PATH}" ]; then
   PATH=/usr/local/; export PATH
else
   PATH=/usr/local/:$PATH; export PATH
fi

if [ -z "${LD_LIBRARY_PATH}" ]; then
   LD_LIBRARY_PATH=/usr/local/:; export LD_LIBRARY_PATH       # Linux
else
   LD_LIBRARY_PATH=/usr/local/::$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
fi

if [ -z "${DYLD_LIBRARY_PATH}" ]; then
   DYLD_LIBRARY_PATH=/usr/local/:; export DYLD_LIBRARY_PATH   # Mac OS X
else
   DYLD_LIBRARY_PATH=/usr/local/::$DYLD_LIBRARY_PATH; export DYLD_LIBRARY_PATH
fi

