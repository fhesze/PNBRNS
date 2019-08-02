function path = getLibPath(OS, Scenario)

mFileName   = mfilename;
path        = which(mFileName);
a           = strfind(path, 'getLibPath.m');
path(a:end) = [];
path        = strcat(path, '/libreac/libbrns_', OS, '/', Scenario, '/');


end