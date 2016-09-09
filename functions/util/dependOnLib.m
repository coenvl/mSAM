function dependOnLib(javalib, url)

if ~exist(javalib, 'file')
    fprintf('Downloading %s for %f', url, javalib)
    websave(javalib, url)
end

javaaddpath(javalib);