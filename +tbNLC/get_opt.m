function val = get_opt(opts, name, default)
%GET_OPT Safe option getter for struct-style options.
    if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
        val = opts.(name);
    else
        val = default;
    end
end

