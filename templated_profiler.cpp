#include<iostream>
#include<cmath>
#include<chrono>
//The ideas can be extended to build a detailed PROFILER...To just profile all the "functions" in a particular C++ application
double variable1;
double variable2;

void initialize(){
    
    for(int i=0; i< 1000000000; ++i){ //Looping for 10^9 times
        variable1 = 2 + i;
    }
}

void another_init(){
    for(int j=0; j<1000000000; ++j){ //Looping for 10^9 times
        variable2 = 2*j;
    }
}
template<typename T>
float time_taken_by_the_function(T(*function)()){
    auto Tstart = std::chrono::high_resolution_clock::now();
    function();
    auto Tend = std::chrono::high_resolution_clock::now();
    auto time_duration = std::chrono::duration_cast<std::chrono::milliseconds>(Tend - Tstart).count();
    std::cout<<"Executing the function took "<< time_duration << "milliseconds" << std::endl;
    return time_duration;
}

int main(){
    std::cout<<"Hi...Check for Profiling...Results coming through"<< std::endl;
    initialize();
    std::cout<<"The initialized value of the variable1 is " << variable1 << std::endl;
    auto time_taken = time_taken_by_the_function(&initialize);
    std::cout<<"Time taken, returned as a RETURN value of the function initalize is " << time_taken << std::endl;
    
        
    another_init();
    std::cout<<"The initialized value of the variable2 is " << variable2 << std::endl;
    auto time_taken2 = time_taken_by_the_function(&another_init);
    std::cout<<"Time taken, returned as a RETURN value of the function another_init is " << time_taken2 << std::endl;

    return 0;
}
