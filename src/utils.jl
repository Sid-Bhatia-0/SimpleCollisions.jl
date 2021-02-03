function Base.findmin(f::Function, itr)
    min_f_value = Inf
    min_element = 0
    min_idx = 0

    for (idx, element) in enumerate(itr)
        f_value = f(element)
        if f_value < min_f_value
            min_f_value = f_value
            min_element = element
            min_idx = idx
        end
    end

    return min_f_value, min_element, min_idx
end

function Base.findmax(f::Function, itr)
    max_f_value = -Inf
    max_element = 0
    max_idx = 0

    for (idx, element) in enumerate(itr)
        f_value = f(element)
        if f_value > min_f_value
            max_f_value = f_value
            max_element = element
            max_idx = idx
        end
    end

    return max_f_value, max_element, max_idx
end

macro pretty_print(name)
    expr = quote
        function Base.show(io::IO, ::MIME"text/plain", object::$(name))
            object_type = typeof(object)
            println("##### ", object_type, " #####")
            for field_name in fieldnames(object_type)
                print(field_name, " => ")
                display(getfield(object, field_name))
            end
        end
    end
end
