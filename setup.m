function setup(opts)

arguments
    opts.surfacefun = 'surfacefun'
    opts.chunkie    = 'chunkie'
end

file = fullfile(opts.surfacefun, 'setup.m');
if ( exist(file, 'file') )
    run(file);
end

file = fullfile(opts.chunkie, 'startup.m');
if ( exist(file, 'file') )
    run(file);
end

end
