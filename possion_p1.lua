function check_array(nodes,str)
    for key,val in ipairs(nodes) do
        print(str, key, table.concat(val, " "))
    end
end 

function check_table(nodes,str)
    for key,val in pairs(nodes) do
        print(str, key, table.concat(val, " "))
    end
end 

function check_matrix(mm, str)
    for row,row_table in ipairs(mm) do
        for col,data in pairs(row_table) do
            print(str,row,col,data)
        end
    end
end 

function check_vector(mm, str)
    for ind,val in ipairs(mm) do
        print(str,ind,val)
    end
end 

function deal_nodes(nodes, Line)
    local n,x,y,z

    n  = io.read("*number")
    for i=1,n do
        n,x,y,z = io.read("*number","*number","*number","*number") 
        if(n) then
            nodes[#nodes+1] = {x,y,z}
            --print(n,x,y,z)
        end
    end
end 

function deal_elements(elements,Line)

    local n,e_type,n_tags,tags,nodes,n_ele,cell

    n  = io.read("*number")
    for i=1,n do
        n_ele= io.read("*number") 
        e_type = io.read("*number")
        n_tags = io.read("*number")
        cell = {n_ele,e_type,n_tags}
        for j=1,n_tags do
            tags = io.read("*number")
            cell[#cell+1]=tags
        end
        if e_type==1 then
            for j=1,2 do
                tags = io.read("*number")
                cell[#cell+1]=tags
            end
        end
        if e_type==2 then
            for j=1,3 do
                tags = io.read("*number")
                cell[#cell+1]=tags
            end
        end 
        elements[n_ele] = cell

    end
end


function method1()

    local nodes         = {}
    local elements      = {}

    local temp = io.input() -- save current file
    io.input("ex0.msh")
    repeat
        local Line = io.read()

        ---[[
        if (Line=="$Nodes") then
            deal_nodes(nodes,Line)
        end
        --]]

        if (Line=="$Elements") then
            deal_elements(elements,Line)
        end

    until Line==nil 

    io.input():close() -- close current file 
    io.input(temp) -- restore previous current file

    check_table(nodes,"nodes: ")
    check_table(elements,"element: ")
end


function read_nodes_2d(nodes, Line)
    local n,x,y,z

    n  = io.read("*number")
    for i=1,n do
        _,x,y,_ = io.read("*number","*number","*number","*number") 

        nodes[i] = {x,y}
    end

end

function read_lines_elements_p1q1(edges, elements, bc_id,Line)

    local n,e_type,n_tags,tag,nodes,n_ele,cell
    local i_edge, i_ele 

    i_edge = 1
    i_ele = 1

    n  = io.read("*number")
    for i=1,n do
        n_ele= io.read("*number") 
        e_type = io.read("*number")
        n_tags = io.read("*number")

        if e_type==1 then

            tag = io.read("*number")
            cell = bc_id[tag] == nil and {} or bc_id[tag]
            cell[#cell+1] = i_edge 
            bc_id[tag] = cell

            for j=2,n_tags do
                --do not save >=2nd  tag right now
                tag = io.read("*number")           
            end

            x,y  = io.read("*number","*numbers")
            cell = {x,y}
            edges[i_edge] = cell

            i_edge = i_edge + 1 -- since edge is read one by one

        elseif e_type==2 then

            for j=1,n_tags do
                tag = io.read("*number")
                --do not save right now
            end

            t1,t2,t3 = io.read("*number","*number","*number")
            cell = {t1,t2,t3}
            elements[i_ele]=cell
            i_ele = i_ele + 1

        elseif e_type==3 then

            for j=1,n_tags do
                tag = io.read("*number")
                --do not save right now
            end

            t1,t2,t3,t4 = io.read("*number","*number","*number","*number")
            cell = {t1,t2,t3,t4}
            elements[i_ele]=cell
            i_ele = i_ele + 1

        end 

    end 
end 

function read_mesh_from_gmsh(filename,nodes,edges,elements,bc_id)

    local temp = io.input() -- save current/default file
    io.input(filename)

    repeat
        local Line = io.read()

        ---[[
        if (Line=="$Nodes") then
            read_nodes_2d(nodes,Line)
        end
        --]]

        ---[[
        if (Line=="$Elements") then
            read_lines_elements_p1q1(edges,elements,bc_id,Line)
        end
        --]]
    until Line==nil 

    io.input():close() -- close current file 
    io.input(temp) -- restore previous current file
    --[[
    check_array(nodes, "node ")
    check_array(edges, "edge ")
    check_array(elements, "elements ")
    check_table(bc_id, "bc_id ")
    --]]
end 

--p1q1 together is because Gmsh in general gives such mesh
function write_p1q1_mesh_vtk(filename, nodes, edges, elements)


    local temp = io.output() -- save current file
    io.output(filename)

    io.write("# vtk DataFile Version 3.1","\n")
    io.write(filename,"\n")
    io.write("ASCII","\n")
    io.write("DATASET UNSTRUCTURED_GRID","\n")
    io.write(string.format("\nPOINTS %5d FLOAT",#nodes),"\n")

    for i=1,(#nodes) do 
        io.write(nodes[i][1],"\t",nodes[i][2],"\t",0,"\n")

    end


    local n_data =0
    for i=1,(#elements) do
        n_data = n_data + 1 + #(elements[i])
    end
    io.write(string.format("\nCELLS %5d %6d",#elements,n_data),"\n")

    for i=1,(#elements) do
        if #(elements[i])==3 then
            io.write(3,"\t",elements[i][1]-1,"\t",elements[i][2]-1,"\t",elements[i][3]-1,"\n")
        elseif #(elements[i])==4 then
            io.write(4,"\t",elements[i][1]-1,"\t",elements[i][2]-1,"\t",elements[i][3]-1,"\t",elements[i][4]-1,"\n")
        end 
    end

    io.write(string.format("\nCELL_TYPES %5d",#elements),"\n")
    for i=1,(#elements) do
        if #(elements[i])==3 then 
            io.write(5,"\t")
        elseif #(elements[i])==4 then 
            io.write(9,"\t")
        end
    end

    io.output():close() -- close current file 
    io.output(temp) -- restore previous current file

end 

function write_p1q1_data_vtk(filename, nodes, edges, elements, u, u_name)


    local temp = io.output() -- save current file
    io.output(filename)

    io.write("# vtk DataFile Version 3.1","\n")
    io.write(filename,"\n")
    io.write("ASCII","\n")
    io.write("DATASET UNSTRUCTURED_GRID","\n")
    io.write(string.format("\nPOINTS %5d FLOAT",#nodes),"\n")

    for i=1,(#nodes) do 
        io.write(nodes[i][1],"\t",nodes[i][2],"\t",0,"\n")

    end

    local n_data =0
    for i=1,(#elements) do
        n_data = n_data + 1 + #(elements[i])
    end
    io.write(string.format("\nCELLS %5d %6d",#elements,n_data),"\n")


    for i=1,(#elements) do
        -- (-1) is because the vtk index of nodes is from 1
        if #(elements[i])==3 then
            io.write(3,"\t",elements[i][1]-1,"\t",elements[i][2]-1,"\t",elements[i][3]-1,"\n")
        elseif #(elements[i])==4 then
            -- note there is only one \n
            io.write(4,"\t",elements[i][1]-1,"\t",elements[i][2]-1,"\t",elements[i][3]-1,"\t",elements[i][4]-1,"\n")
        end 
    end

    io.write(string.format("\nCELL_TYPES %5d",#elements),"\n")
    for i=1,(#elements) do
        if #(elements[i])==3 then 
            io.write(5,"\t")
        elseif #(elements[i])==4 then 
            io.write(9,"\t")
        end
    end

    io.write(string.format("POINT_DATA  %5d",#u),"\n")
    io.write(string.format("SCALARS %s FLOAT",u_name),"\n")
    io.write("LOOKUP_TABLE default\n")
    for i=1,(#u) do 
        io.write(u[i],"\n")

    end


    io.output():close() -- close current file 
    io.output(temp) -- restore previous current file

end 

function initialize_matrix(mesh,mm)
    local ele = {}
    local row_table = {}
    for i=1,#mesh.elements do
        ele = mesh.elements[i]
        for j=1,#ele do
            row_table = mm[ele[j]]==nil and {} or mm[ele[j]]
            for k=1,#ele do 
                --number of elements at a node or edge
                --row_table[ele[k]] = row_table[ele[k]]==nil and 1 or (row_table[ele[k]]+1)
                row_table[ele[k]] = 0
                --print(ele[j],ele[k],row_table[ele[k]])
            end
            mm[ele[j]] = row_table --since nil is not {}
        end
    end 
end 

function get_triangle_J(nodes)
    local J,invJ,detJ
    J={{},{}}
    invJ={{},{}}
    J[1][1]=nodes[2][1]-nodes[1][1]
    J[2][1]=nodes[2][2]-nodes[1][2]
    J[1][2]=nodes[3][1]-nodes[1][1]
    J[2][2]=nodes[3][2]-nodes[1][2]
    detJ = J[1][1]*J[2][2]-J[1][2]*J[2][1]
    invJ[1][1]=J[2][2]/detJ
    invJ[1][2]=-J[1][2]/detJ
    invJ[2][1]=-J[2][1]/detJ
    invJ[2][2]=J[2][2]/detJ
    return J,invJ,detJ
end 

function add_local_matrix_to_global(local_MC,dofs,MC)
    for i=1,#dofs do
        local indx_i = dofs[i]
        for j=1,#dofs do
            local indx_j = dofs[j]
            MC[indx_i][indx_j] = MC[indx_i][indx_j] + local_MC[i][j]
        end
    end
end 

function get_quadrature_triangle(method)
    local quad,np

    --if method=="CENTROID" then

    quad={

        points = { {0.33333333333333333333,  0.33333333333333333333} },
        weight = {1}
    }

    np = 1

    if method=="GAUSS4X4" then

        quad={
            points = { 
                {0.0571041961, 0.06546699455602246},
                {0.2768430136, 0.05021012321401679},
                {0.5835904324, 0.02891208422223085},
                {0.8602401357, 0.009703785123906346},
                {0.0571041961, 0.3111645522491480},
                {0.2768430136, 0.2386486597440242},
                {0.5835904324, 0.1374191041243166},
                {0.8602401357, 0.04612207989200404},
                {0.0571041961, 0.6317312516508520},
                {0.2768430136, 0.4845083266559759},
                {0.5835904324, 0.2789904634756834},
                {0.8602401357, 0.09363778440799593},
                {0.0571041961, 0.8774288093439775},
                {0.2768430136, 0.6729468631859832},
                {0.5835904324, 0.3874974833777692},
                {0.8602401357, 0.1300560791760936}
            },


            weight = { 
                0.04713673637581137,
                0.07077613579259895,
                0.04516809856187617,
                0.01084645180365496,
                0.08837017702418863,
                0.1326884322074010,
                0.08467944903812383,
                0.02033451909634504,
                0.08837017702418863,
                0.1326884322074010,
                0.08467944903812383,
                0.02033451909634504,
                0.04713673637581137,
                0.07077613579259895,
                0.04516809856187617,
                0.01084645180365496
            }    
        }

        np=16
    elseif method=="STRANG1" then 
        quad =
        {
            points = { 
                {0.66666666666666666667,  0.16666666666666666667},
                {0.16666666666666666667,  0.66666666666666666667},
                {0.16666666666666666667,  0.16666666666666666667}
            },

            weight = { 0.33333333333333333333,
            0.33333333333333333333,
            0.33333333333333333333
        }

    }

    np = 3

elseif method == "STRANG3" then 

    quad = {
        points= { 
            {0.33333333333333333333, 0.33333333333333333333},
            {0.60000000000000000000, 0.20000000000000000000},
            {0.20000000000000000000, 0.60000000000000000000},
            {0.20000000000000000000, 0.20000000000000000000}
        },

        weight = {  -0.56250000000000000000,
        0.52083333333333333333,
        0.52083333333333333333,
        0.52083333333333333333
    }

}

np = 4

 end
return quad,np
end 


function assemble_mass_matrix(mesh,MC)
    local local_mc = {}
    local dofs,nodes,detJ,invJ,J

    for i=1,#mesh.elements do
        dofs = mesh.elements[i]
        nodes = {}
        for j=1,#dofs do 
            nodes[#nodes+1] = mesh.nodes[dofs[j]]  
        end 

        J,invJ,detJ = get_triangle_J(nodes)

        local_mc = {{},{},{}}
        for j=1,#dofs do
            for k=1,#dofs do
                local_mc[j][k] = j==k and detJ/12 or   detJ/24
            end
        end
        add_local_matrix_to_global(local_mc,dofs,MC)
    end 
end


function add_local_vector_to_global(local_vec,dofs,vec)
    for i=1,#dofs do
        local indx_i = dofs[i]
        vec[indx_i] = vec[indx_i] + local_vec[i]
    end
end 

function assemble_lumped_mass_matrix(mesh,ML)
    local local_ml = {}
    local dofs,nodes,detJ,J,invJ

    for i=1,#mesh.elements do
        dofs = mesh.elements[i]
        nodes = {}
        for j=1,#dofs do 
            nodes[#nodes+1] = mesh.nodes[dofs[j]]  
        end 

        J,invJ,detJ = get_triangle_J(nodes)
        for j=1,#dofs do
            local_ml[j] = detJ/6 
        end
        add_local_vector_to_global(local_ml,dofs,ML)
    end     
end

--the product of (n1 x n2) matrix mm1 and (n2 x n3) matrixmm2
function matmul(mm1,mm2,n1,n2,n3)
    local mmres={}
    for i=1,n1 do
        mmres[i]={}
        for j=1,n3 do
            for k=1,n2 do 
                mmres[i][j] =  mmres[i][j]==nil and mm1[i][k]*mm2[k][j] or mmres[i][j]+mm1[i][k]*mm2[k][j]
            end
        end
    end
    return mmres
end 

function sparse_vmult(mm,vec)
    local vecres={}
    local row={}
    
    for i=1,#mm do
        vecres[i]=0
        row = mm[i]
        for col,val in pairs(row) do
            vecres[i] = vecres[i] + val*vec[col]
        end
    end

    return vecres
end

function vecdot(vec1,vec2)
    local res
    res = 0
    for i=1,#vec1 do
        res = res + vec1[i]*vec2[i]
    end
    return res
end 

function set_spmat(mm,cc)
    for i=1,#mm do
        for j,_ in pairs(mm[i]) do 
            mm[i][j] = cc
        end
    end 
end 

function assemble_stiffness_matrix(mesh,ss)
    local local_ss = {}
    local dofs,nodes,detJ,invJ,J
    local shape_grad
    
    set_spmat(ss,0)
     
    for i=1,#mesh.elements do
        dofs = mesh.elements[i]
        nodes = {}
        for j=1,#dofs do 
            nodes[j] = mesh.nodes[dofs[j]]  
        end 
    
        J,invJ,detJ = get_triangle_J(nodes)
        shape_grad = matmul({{-1,-1},{1,0},{0,1}},invJ,3,2,2)
           
        local_ss = {}
        for j=1,#dofs do
            local_ss[j]={}
            for k=1,#dofs do
                local_ss[j][k] = vecdot(shape_grad[j],shape_grad[k])*detJ*0.5 --const*area
            end
        end
        
        add_local_matrix_to_global(local_ss,dofs,ss)
    end 
end

function get_p1_fe_base(nodes,quad_point,weight)
    local shape_grad,shape_value,JxW
    local J,invJ,detJ,cell_point
    J,invJ,detJ = get_triangle_J(nodes)    
    shape_grad = matmul({{-1,-1},{1,0},{0,1}},invJ,3,2,2)
    shape_value= {1-quad_point[1]-quad_point[2],quad_point[1],quad_point[2]}
    JxW = 0.5*detJ*weight
    cell_point={}
    cell_point[1]=J[1][1]*quad_point[1]+J[1][2]*quad_point[2]+nodes[1][1]
    cell_point[2]=J[2][1]*quad_point[1]+J[2][2]*quad_point[2]+nodes[1][2]
    return shape_value,shape_grad,JxW,cell_point
end 

function assemble_stiffness_matrix_M2(mesh,ss)
    local local_ss = {}
    local dofs,nodes,detJ,invJ,J
    local shape_grad,shape_value,JxW
    
    local quad,np
    quad,np = get_quadrature_triangle("STRANG3")
    
    set_spmat(ss,0)
     
    for i=1,#mesh.elements do
        dofs = mesh.elements[i]
        nodes = {}
        for j=1,#dofs do 
            nodes[j] = mesh.nodes[dofs[j]]  
        end 
    
        local_ss = {{},{},{}}
        for ip=1,np do 
            shape_value,shape_grad,JxW = get_p1_fe_base(nodes,quad.points[ip],quad.weight[ip])
                        
            for j=1,#dofs do
                for k=1,#dofs do
                    local_ss[j][k] = local_ss[j][k]==nil and vecdot(shape_grad[j],shape_grad[k])*JxW or 
                                     local_ss[j][k]+ vecdot(shape_grad[j],shape_grad[k])*JxW
                end
            end
            
        end
        add_local_matrix_to_global(local_ss,dofs,ss)
    end 
end

function assemble_mass_matrix_M2(mesh,mc)
    local local_mc = {}
    local dofs,nodes,detJ,invJ,J
    local shape_grad,shape_value,JxW
    
    local quad,np
    quad,np = get_quadrature_triangle("STRANG3")
    
    set_spmat(mc,0)
     
    for i=1,#mesh.elements do
        dofs = mesh.elements[i]
        nodes = {}
        for j=1,#dofs do 
            nodes[j] = mesh.nodes[dofs[j]]  
        end 
    
        local_mc = {{},{},{}}
        for ip=1,np do 
            shape_value,shape_grad,JxW = get_p1_fe_base(nodes,quad.points[ip],quad.weight[ip])
                        
            for j=1,#dofs do
                for k=1,#dofs do
                    local_mc[j][k] = local_mc[j][k]==nil and shape_value[j]*shape_value[k]*JxW or 
                                     local_mc[j][k]+shape_value[j]*shape_value[k]*JxW
                end
            end
            
        end
        add_local_matrix_to_global(local_mc,dofs,mc)
    end 
end

function set_vector(vec,cc)
    for i=1,#vec do
        vec[i]=cc
    end
end 

function assemble_lumped_mass_matrix_M2(mesh,ml)
    local local_ml = {}
    local dofs,nodes,detJ,invJ,J
    local shape_grad,shape_value,JxW
    
    local quad,np
    quad,np = get_quadrature_triangle("STRANG1")
    
    set_vector(ml,0)
     
    for i=1,#mesh.elements do
        dofs = mesh.elements[i]
        nodes = {}
        for j=1,#dofs do 
            nodes[j] = mesh.nodes[dofs[j]]  
        end 

        local_ml = {}
        for ip=1,np do 
            shape_value,shape_grad,JxW = get_p1_fe_base(nodes,quad.points[ip],quad.weight[ip])

            for j=1,#dofs do
                local_ml[j] = local_ml[j]==nil and shape_value[j]*JxW or 
                local_ml[j]+shape_value[j]*JxW

            end

        end
        add_local_vector_to_global(local_ml,dofs,ml)
    end 
end


function func_f(pt)
    local x=pt[1]
    local y=pt[2]
    return 32*y*(1-y)+32*x*(1-x)
end

function func_solution(pt)
    local x=pt[1]
    local y=pt[2]
    return 16*x*(1-x)*y*(1-y)
end 

function assemble_rhs_by_function(mesh,func,rhs)
    local local_rhs = {}
    local dofs,nodes,detJ,invJ,J
    local shape_grad,shape_value,JxW
    local f_at_p,cell_point
    local quad,np
    quad,np = get_quadrature_triangle("STRANG1")

    set_vector(rhs,0)

    for i=1,#mesh.elements do
        dofs = mesh.elements[i]
        nodes = {}
        for j=1,#dofs do 
            nodes[j] = mesh.nodes[dofs[j]]  
        end 

        local_rhs = {}
        for ip=1,np do 
            shape_value,shape_grad,JxW,cell_point = get_p1_fe_base(nodes,quad.points[ip],quad.weight[ip])
            f_at_p = func(cell_point)
            for j=1,#dofs do
                local_rhs[j] = local_rhs[j]==nil and f_at_p*shape_value[j]*JxW or 
                local_rhs[j]+f_at_p*shape_value[j]*JxW
            end

        end
        add_local_vector_to_global(local_rhs,dofs,rhs)
    end 
end

function test1()
    --global variables

    --local variables
    local mesh = {}
    mesh.nodes = {}
    mesh.edges = {}
    mesh.elements = {}
    mesh.bc_id = {}--begin with 0

    --start time    
    start_time = os.clock()

    --create a mesh
    read_mesh_from_gmsh("square2x2.msh",mesh.nodes,mesh.edges,mesh.elements,mesh.bc_id)
    write_p1q1_mesh_vtk("mesh.vtk",mesh.nodes,mesh.edges,mesh.elements)


    --create a vector 
    local u={}
    for i=1,#mesh.nodes do
        --u[i]= mesh.nodes[i][1]*(1.0-mesh.nodes[i][1])*mesh.nodes[i][2]*(1.0-mesh.nodes[i][2])
        u[i] = 0
    end 

    write_p1q1_data_vtk("data.vtk",mesh.nodes,mesh.edges,mesh.elements,u,"u")

    --assemble the consistent mass matrix
    local MC = {}
    initialize_matrix(mesh,MC)
    assemble_mass_matrix(mesh,MC)
    check_matrix(MC,"mass matrix: ")
    assemble_mass_matrix_M2(mesh,MC)
    check_matrix(MC,"mass matrix-M2: ")

    assemble_lumped_mass_matrix(mesh,u)
    write_p1q1_data_vtk("ML.vtk",mesh.nodes,mesh.edges,mesh.elements,u,"ml")
    check_vector(u,"ml: ")
    assemble_lumped_mass_matrix_M2(mesh,u)
    check_vector(u,"ml-M2: ")
    assemble_stiffness_matrix(mesh,MC)
    check_matrix(MC,"stiff matrix: ")
    assemble_stiffness_matrix_M2(mesh,MC)
    check_matrix(MC,"stiff matrix-M2: ")
    
    --assemble right hand side
    local rhs={}
    for i=0,#mesh.nodes do
        rhs[i]=0
    end
    assemble_rhs_by_function(mesh,func_f,rhs)
    write_p1q1_data_vtk("rhs.vtk",mesh.nodes,mesh.edges,mesh.elements,rhs,"rhs")
    
    --check time
    end_time = os.clock()
    elapsed_time = end_time-start_time
    print(elapsed_time)

    return elapsed_time
end 

function func_tag_to_value(tag)
    local value 
    value = 0
    if(tag == 5)then
        value = 1
    end
    
    return value
end

function eliminate_row(row_index,value,MM,rhs)
    local row = MM[row_index]
    for col_index,_ in pairs(row) do 
        if(row_index==col_index) then
            row[col_index] = 1
            rhs[row_index] = value
        else
            rhs[col_index]=rhs[col_index]-MM[col_index][row_index]*value
            MM[col_index][row_index] = 0
            row[col_index] = 0
        end
    end 
end

function apply_bc(mesh,SS,func,rhs)
    --avoid eliminate twice for each nodes since 
    --bc_id is labeled for edges
    local visited={}
    --since bc_id begins with 0, one has to use pairs instead of ipairs
    for tag,lines in pairs(mesh.bc_id) do 
        if(tag ~=5) then--natural boundary condition
            for i=1,#lines do 
                --print(tag,i,lines[i],mesh.edges[lines[i]][1],mesh.edges[lines[i]][2])
                if(visited[mesh.edges[lines[i]][1]]==nil) then 
                    eliminate_row(mesh.edges[lines[i]][1],func_tag_to_value(tag),SS,rhs)
                    visited[mesh.edges[lines[i]][1]] = 1    
                end
                if(visited[mesh.edges[lines[i]][2]]==nil) then 
                    eliminate_row(mesh.edges[lines[i]][2],func_tag_to_value(tag),SS,rhs)
                    visited[mesh.edges[lines[i]][2]] = 1    
                end
            end
        end
    end
    
end

function initialize_vector(nn,cc)
    local u={}
    for i=1,nn do
        u[i] = cc
    end
    return u
end 

function assign_vector(xx,bb)
    for ind,_ in pairs (xx) do
        xx[ind] = bb[ind]
    end 
end 

function vector_sadd(a,u,b,v)
    for ind,_ in pairs (u) do
        u[ind] = a*u[ind]+b*v[ind]
    end 
end

--xx has initial guess x0
function solve_Axb_cg(AA,xx,bb,epsilon,max_it)
    local res = sparse_vmult(AA,xx)
    vector_sadd(-1,res,1,bb)
    local rho0 = vecdot(res,res)
    local k=1
    local p=initialize_vector(#xx,0)
    local beta,w,rho1,rho2,a
    rho1 = rho0
    rho2 = rho0
    repeat 
        if(k==1)then
            assign_vector(p,res)
        else
            beta = rho2/rho1
            vector_sadd(beta,p,1,res)
        end
        w = sparse_vmult(AA,p)
        a = rho2/vecdot(p,w)
        vector_sadd(1,xx,a,p)
        vector_sadd(1,res,-a,w)

        rho1 = rho2
        rho2 = vecdot(res,res)
        --print (rho2,rho1,rho0)
        k = k+1
    until (rho2<rho0*epsilon or k>max_it)
end

function get_L2_error_between_vec_and_func(mesh,uu,ff)
       
    local dofs,nodes,detJ,invJ,J
    local shape_grad,shape_value,JxW
    local f_at_p,u_at_p,cell_point
    local quad,np
    quad,np = get_quadrature_triangle("")

    local err=0
    
    for i=1,#mesh.elements do
        dofs = mesh.elements[i]
        nodes = {}
        for j=1,#dofs do 
            nodes[j] = mesh.nodes[dofs[j]]  
        end 

        for ip=1,np do 
            shape_value,shape_grad,JxW,cell_point = get_p1_fe_base(nodes,quad.points[ip],quad.weight[ip])
            --serious error: use point in real cell
            f_at_p = ff(cell_point)
            
            u_at_p = 0
            for j=1,#dofs do
                u_at_p = u_at_p + uu[dofs[j]]*shape_value[j]
            end
            
            err = err + (f_at_p-u_at_p)*(f_at_p-u_at_p)*JxW
        end
    end 

    return math.sqrt(err)
end
function test2()
        
    --local variables
    local mesh = {}
    mesh.nodes = {}
    mesh.edges = {}
    mesh.elements = {}
    mesh.bc_id = {}--begin with 0

    --start time
    start_time = os.clock()

    --create a mesh
    read_mesh_from_gmsh("square512x512.msh",mesh.nodes,mesh.edges,mesh.elements,mesh.bc_id)
    write_p1q1_mesh_vtk("mesh.vtk",mesh.nodes,mesh.edges,mesh.elements)

        
    --initialize rhs
    local rhs={}
    for i=1,#mesh.nodes do
        rhs[i]=func_solution(mesh.nodes[i])
        --rhs[i] = 0
    end
    write_p1q1_data_vtk("true.vtk",mesh.nodes,mesh.edges,mesh.elements,rhs,"uu")

    --interpolation error    
    local err = get_L2_error_between_vec_and_func(mesh,rhs,func_solution)
    print(#mesh.nodes,err)
    
    --assemble right hand side
    assemble_rhs_by_function(mesh,func_f,rhs)
    write_p1q1_data_vtk("rhs.vtk",mesh.nodes,mesh.edges,mesh.elements,rhs,"rhs")
     
    --assemble stiffness matrix
    local SS={}
    initialize_matrix(mesh,SS)
    assemble_stiffness_matrix_M2(mesh,SS)
        
    --apply essential boundary condition
    --use mesh.bc_id
    --check_matrix(SS,"stiff matrix before applying bc: ")    
    apply_bc(mesh,SS,func_tag_to_value,rhs)
    --check_matrix(SS,"stiff matrix after applying bc: ")

    --initialize solution
    local u={}
    for i=1,#mesh.nodes do
        u[i] = 0
    end
    --solve linear system using CG
    solve_Axb_cg(SS,u,rhs,1e-6,1000)
    
    write_p1q1_data_vtk("uu.vtk",mesh.nodes,mesh.edges,mesh.elements,u,"u")
    
    --compute the error
    err = get_L2_error_between_vec_and_func(mesh,u,func_solution)
    print(#mesh.nodes,err)

    --check time
    end_time = os.clock()
    elapsed_time = end_time-start_time
    print(elapsed_time)

    return elapsed_time

end

test2(...)


