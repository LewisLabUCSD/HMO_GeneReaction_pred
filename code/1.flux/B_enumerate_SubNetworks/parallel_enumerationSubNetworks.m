%Here below two examples on how we generated all the subnetworks

solution={};
parfor z=2485:2491
      [sol]=enumerate_AllsubNetworks(z);
      solution{z}=sol;
end
save combi_2485-2491

solution={};
parfor z=2786:2793
      [sol]=enumerate_AllsubNetworks(z);
      solution{z}=sol;
end
save combi_2786-2793





