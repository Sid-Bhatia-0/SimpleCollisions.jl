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
