%Finite element solver for one dimensional poroelasticity
%Adilan W Mahdiyasa
%First version 2024

%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% Input Data Main Program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
    input_data               = input_1D;     
    time_max                 = 1000;                       %Maximum simulation time 
    dt                       = 0.1;                        %Time interval
    initial_condpp           = 1e5;                        %Initial condition of pore water pressure
    number_nodes             = input_data.control(1);      %Total nodes used in the simulation
    number_element           = input_data.control(2);      %Total element used in the simulation
    nodaldof                 = input_data.control(3);      %Degree of freedom
    nodes_perelement         = input_data.control(4);      %Total nodes per element
    number_material          = input_data.control(5);      %Total types of material
    numeq                    = number_nodes*nodaldof;
    totalelementdof          = nodaldof*nodes_perelement;
%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% Initial Array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
    KG          = zeros(numeq,numeq);   %Global Stiffness Matrix
    CG          = zeros(numeq,numeq);   %Global Damping Matrix
    F           = zeros(numeq,1);       %Global loads or force
    solution    = zeros(numeq,1);       %Solution of poroelasticity
    
%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% Main Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
% Assemble Global Stiffness KG and Damping Matrix CG 
    for iel=1:number_element;  
         [local_data]        = localize_1D(iel,input_data);
         Element_No          = iel;
         Material_Array      = local_data.mater;
         Coordinate_Array    = local_data.coords;
         dofAddress_Array    = local_data.dofs;
         Material_Type       = input_data.EL(iel,2);
         Element_Type        = input_data.mater(Material_Type,1);
         [Kl,Cl]             = element_1D( local_data );

 %        Creating local stiffness matrix 
         Element_No         = iel;
         LocalStiffness     = Kl;
         LocalDamping       = Cl;

 %        Assembling Global Stiffness and Global Damping matrix 
         [KG]               = assemble_1D (KG, Kl, local_data, totalelementdof); %Assemble global stiffness matrix KG from local stiffness matrix Kl
         [CG]               = assemble_1D (CG, Cl, local_data, totalelementdof); %Assemble global damping matrix GG from local damping matrix Cl
         
         Element_No         = iel;
         Global_Stiffness   = KG;
         Global_Damping     = CG;
            
    end
%         Assembling load and boundary conditions
        [F]                 = vectorstretch_1D (input_data.LOAD, number_nodes, 1+nodaldof);  %Creating vector column of load
        [con]               = vectorstretch_1D (input_data.CON,  number_nodes, 1+nodaldof) ;  %Creating vector column of boundary condition
        solution            = F.* (con); 
        F                   = F.* (1-con); 
        F                   = F - KG*solution; 

%         Creating system of equations                           
        kg                      = KG;                       %Temporary matrix for global stiffness
        cg                      = CG;                       %Temporary matrix for global damping
        f                       = F;                        %Temporary matrix for global load
        sol_temp                = zeros(numeq,1);
        sol_temp(2:2:numeq)     = 1;
        sol_temp                = sol_temp*initial_condpp;


    for ieq = numeq: -1: 1;
        if (con(ieq,1) == 1)
            kg(:,ieq)           =[]; 
            kg(ieq,:)           =[]; 
            cg(:,ieq)           =[]; 
            cg(ieq,:)           =[];                          
            f(ieq,:)            =[]; 
            sol_temp(ieq,:)     =[];
        end
    end

%        Solving system of equations
    for itime = 1:time_max;                   
                sol_temp        = (kg + (1/dt)*cg)\(f + (1/dt)*(cg*sol_temp));
                f               = f*0;
    end
icount=0;                                                 
    for ieq = 1:numeq;
        if (con(ieq,1) == 0) 
            icount          = icount + 1; 
            solution(ieq,1) = sol_temp(icount,1);
        end
    end
    solution;                                                       
for i=1:number_nodes 
        displacement(i)     =solution((2*i)-1);     %Solution of displacement
        pore_pressure(i)    =solution(2*i);         %Solution of pore pressure
end


 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  
    '***************** Thank You ******************* '