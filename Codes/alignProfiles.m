function [alignedProfile1 alignedProfile2 shortlength]=alignProfiles(LC_profile1,LC_profile2)

  if (nnz(LC_profile1)>1)
 Startpoint = 1;
 while LC_profile1(Startpoint,1)==0
    Startpoint = Startpoint + 1;
 end
 Endpoint=length(LC_profile1);
 while LC_profile1(Endpoint,1)==0
    Endpoint = Endpoint - 1;
 end
 LC_profile1=LC_profile1(Startpoint:Endpoint,1);
 end

 if (nnz(LC_profile2)>1)
 Startpoint = 1;
 while LC_profile2(Startpoint,1)==0
    Startpoint = Startpoint + 1;
 end
 Endpoint=length(LC_profile2);
 while LC_profile2(Endpoint,1)==0
    Endpoint = Endpoint - 1;
 end
 LC_profile2=LC_profile2(Startpoint:Endpoint,1);
 end
 
 
  l1=nnz(LC_profile1);
  l2=nnz(LC_profile2);
  if (l2>=l1/2 && l1>1)
      if l1>l2
      shorterProfile=LC_profile2;
      longerProfile=LC_profile1;
      shortlength=l2;
      longlength=l1;
  else
       shorterProfile=LC_profile1;
       longerProfile=LC_profile2;
       shortlength=l1;
       longlength=l2;
  end     
  longerProfile=[longerProfile ;zeros(shortlength,1)];
  
  flipedProfile=shorterProfile(shortlength:-1:1);
  
  matchedFilterOutput=conv(longerProfile,flipedProfile,'full');
 
  [~, maxposi]=max(matchedFilterOutput);
  
  if maxposi> longlength
      longEnd=longlength;
     
  else
      longEnd=maxposi;
  end
   endZeros=maxposi-longlength;
  
  if maxposi-shortlength+1<1
      longStart=1;
  else
      longStart= maxposi-shortlength+1;
  end
  startZeros=1-(maxposi-shortlength+1);
  
  alignedLongerProfile=longerProfile(longStart:longEnd);
  if endZeros>0
      alignedLongerProfile=[alignedLongerProfile; zeros(endZeros,1)];
  end
  if startZeros>0
      alignedLongerProfile=[zeros(startZeros,1);alignedLongerProfile; ];
  end
  if l1>l2
      alignedProfile1=alignedLongerProfile;
      alignedProfile2=shorterProfile;
  else
      alignedProfile2=alignedLongerProfile;
      alignedProfile1=shorterProfile;
  end  
  else
   alignedProfile1=LC_profile1;
   alignedProfile2=LC_profile2;  
   shortlength=1;
  end 

  
  