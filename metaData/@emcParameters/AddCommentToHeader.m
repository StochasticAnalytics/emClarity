function [ ] = AddCommentToHeader(obj,comment_to_add)
%   Detailed explanation goes here

  % 
% 	if (comment_to_add.StartsWith("#") == false)
  if isempty(regexp(comment_to_add,'^#*','once'))
	
		comment_to_add = "# " + comment_to_add;
    
  end

	comment_to_add = strtrim(comment_to_add);
  
	obj.header_comments.Add(comment_to_add);
  
end

