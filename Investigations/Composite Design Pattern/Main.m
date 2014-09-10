function Main
    containers = OrchestraCompositeCl();
    subContainers = OrchestraCompositeCl();
    obj = OrchestraCl(1);
    subContainers = subContainers.add(obj);
    subContainers = subContainers.add(OrchestraCl(2));
    subContainers = subContainers.add(OrchestraCl(3));
    containers = containers.add(OrchestraCl(4));
    containers = containers.add(subContainers);
    containers = containers.add(containers);
    
%     subContainers.Calculate()
    result = containers.Calculate()
end
