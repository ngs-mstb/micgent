{
  data : {
    local root = std.extVar("data"),
    cge : {
        // Defining identical variable names in enclosing scopes
        // leads to infinite recursion, and there is no way to do relative 
        // dereferencing from "one scope up". Therefore, we use
        // variables with unique names
        local root_cge = root + "/cge",
        plasmid_db : {
            local root_plasmid_db = root_cge + "/plasmid_finder/2016-01-04/plasmidfinder",
            seq : root_plasmid_db + "/plasmid_database.fsa"
        },
        res_db : {
            local root_res_db = root_cge + "/res_finder/2016-03-01/database",
            // in some versions there is an 'all.fsa' file or similar -
            // need to exclude it
            seq : root_res_db + "/*.fsa"
        }
    }
  },
  wrapper : null
}
