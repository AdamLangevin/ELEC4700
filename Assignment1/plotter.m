function plotter(Limits)

global Vx Vy x y MarkerSize T time Elect

if isempty(Limits)
    Limits = Limits;
end
hold on
plot(Elect(:,1), Elect(:,2),'markers',MarkerSize);
axis(Limits);

end

