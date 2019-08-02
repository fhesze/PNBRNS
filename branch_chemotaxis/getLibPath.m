function path = getLibPath(OS)

mFileName   = mfilename;
path        = which(mFileName);
a           = strfind(path, 'getLibPath.m');
path(a:end) = [];
path        = strcat(path, 'libreac/libbrns_', OS, '/');


end