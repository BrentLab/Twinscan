#!/bin/sh


if [ -x /usr/bin/pkg-config ] || [ -x /usr/local/bin/pkg-config ]; then
	echo -n "GLIB_CFLAGS="
    echo -n `pkg-config --cflags glib-2.0`
    echo " -DHAS_GLIB"

	echo -n "GLIB_LFLAGS="
    echo -n `pkg-config --libs glib-2.0`
    echo " -static -DHAS_GLIB"
else
	echo "GLIB_CFLAGS="
	echo "GLIB_LFLAGS="
fi
