function grade=getGrade(score)

if score >= 89
    grade = 'Excellent';
elseif score >= 77
    grade = 'Great';
elseif     score >= 65
    grade = 'Satisfactory';
else
    grade = 'Try Again';
 end
    
end
