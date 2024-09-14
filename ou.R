rm(list = ls()) # cleans environment 

simulate_ou_step = function(x,y,alpha,sigma,mu)
{
  
  new_x = x + alpha*(mu - x) + rnorm(1,0,sigma)
  new_y = y + alpha*(mu - y) + rnorm(1,0,sigma)
  return(list(new_x,new_y))
}

multiple_simulate_ou_step = function(x,y,alpha,sigma,mu)
{
  new_x = x
  new_y = y
  for (i in 1:length(x)) # for each turtles x and y position
  {
    new_position = simulate_ou_step(x[i],y[i],alpha,sigma,mu)
    new_x[i] = new_position[[1]]
    new_y[i] = new_position[[2]]
  }
  return(list(new_x,new_y))
}

simulate_turtles = function(alpha,sigma,mu,T,n)
{
  positions = matrix(0,(T+1)*n,2) 
  for (i in 1:n) # for each turtle , set their initial position to either 0 or 1
  {
    positions[i,1] = round(runif(1)) 
    positions[i,2] = round(runif(1)) 
  }
  
  for (t in 1:T) # for each time-step 
  { 
    new_positions = multiple_simulate_ou_step(positions[(t-1)*n+(1:n),1],positions[(t-1)*n+(1:n),2],alpha,sigma,mu)
    positions[t*n+(1:n),1] = new_positions[[1]] 
    positions[t*n+(1:n),2] = new_positions[[2]]
  }
  
  return(positions)
}

final_positions = simulate_turtles(0.03, 1.22, 50, 100, 10)




