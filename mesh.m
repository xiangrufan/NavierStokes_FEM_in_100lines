
classdef mesh <handle
    %BOUNDARY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Position;
        TriangleIndex;
        Edge;
        Boundaries;
    end
    
    methods % stores initiator function 
        function MESH=mesh(p,t)
            MESH.Position=p;
            MESH.TriangleIndex=t; 
            MESH.Edge=boundedges(p,t);
        end
        
        
        
    end
    methods % stores BC functions
        
        function bc = rec_selector(MESH,LowerLeft,UpperRight)
            p=MESH.Position;
            t=MESH.TriangleIndex;
            e=unique(MESH.Edge);
            [num_edge,~]=size(e);
            i=1;
%             bc=zeros(1,2);
            for ele=1:num_edge
          EdgPos=p(e(ele),:);
              
            if (  EdgPos(1)>=LowerLeft(1) && EdgPos(2)>=LowerLeft(2) &&  EdgPos(1)<=UpperRight(1) && EdgPos(2)<=UpperRight(2))
%  check point here               bc(i,:)=EdgPos ;
bc(i)=e(ele);
                i=i+1;
            end
            
            end
            
        end
     end
    
end

