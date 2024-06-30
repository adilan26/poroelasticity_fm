function [local_data] = localize_1D (iel, in_data)

% Localize global data at the local element level
% Define global system variables

  number_nodes           = in_data.control(1);  %Total nodes
  number_element         = in_data.control(2);  %Total element
  nodaldof               = in_data.control(3);  %Degree of freedom
  nodes_perelement       = in_data.control(4);  %Total nodes per element
  number_material        = in_data.control(5);  %Material type

  local_data.mater   = zeros(1,5);
  local_data.coords  = zeros(nodes_perelement,1);
  local_data.dofs    = zeros(1,nodes_perelement*nodaldof);
  material_number    = in_data.EL(iel,2);  
  local_data.mater   = in_data.mater(material_number,1:4);                      
        for inodes=1:nodes_perelement                                              
            local_data.coords(inodes,1:1) = in_data.ND(in_data.EL(iel,2+inodes),2:1+1) ;
        end
        
        for inodes=1:nodes_perelement
            globalnode = in_data.EL(iel,2+inodes);
            startdof   = (globalnode-1)*nodaldof ;                              
            for idof=1:nodaldof
            local_data.dofs(1,(inodes-1)*nodaldof+idof) = startdof+idof;  
            end
        end