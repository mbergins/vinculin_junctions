function mkdir_no_err(dir)

if not(exist(dir,'dir'))
    mkdir(dir);
end