#include <dynsyslab/workflow/Queue.hpp>

#include <chrono>
#include <filesystem>
#include <future>
#include <iostream>
#include <stdexcept>
#include <thread>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

dynsyslab::workflow::json::Value solve(double stop, double step) {
    using dynsyslab::workflow::json::parse;
    auto job = parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"queue-test",
      "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
      "time":{"start":0,"stop":1,"output_step":0.01,"coordinate_system":"automatic"},
      "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"seed":0,
        "values":[[-1,0.1]],"ranges":[]},
      "controller":{"kind":"none","settings":{},
        "target":{"kind":"constant","state":[-1,0]}},
      "solver":{"method":"RK4"},
      "archive":{"enabled":true,"data_role":"free_training"}
    })");
    job.at("time").at("stop") = stop;
    job.at("time").at("output_step") = step;
    return job;
}

}  // namespace

int main() {
    const auto root=std::filesystem::temp_directory_path()/"dynsyslab-workflow-queue-test";
    std::error_code ignored;std::filesystem::remove_all(root,ignored);
    try {
        dynsyslab::workflow::JobQueue queue(root);
        const auto completedId=queue.submit(solve(0.02,0.01));
        const auto completed=queue.runNext("queue-test");
        require(completed.state==dynsyslab::workflow::QueueRunResult::State::completed,
                "short queued solve did not complete");
        require(completed.job&&completed.job->id==completedId&&completed.run,
                "completed queued result lost identity or workflow result");
        require(!completed.run->archives.empty(),"completed queued solve has no archive");

        const auto cancelledId=queue.submit(solve(1000.0,0.001));
        auto running=std::async(std::launch::async,[&]{return queue.runNext("queue-cancel-test");});
        const auto deadline=std::chrono::steady_clock::now()+std::chrono::seconds(5);
        while(std::chrono::steady_clock::now()<deadline){
            const auto record=queue.find(cancelledId);
            if(record&&record->status==dynsyslab::core::JobStatus::running)break;
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
        (void)queue.cancel(cancelledId);
        require(running.wait_for(std::chrono::seconds(10))==std::future_status::ready,
                "queued cancellation was not observed promptly");
        const auto cancelled=running.get();
        require(cancelled.state==dynsyslab::workflow::QueueRunResult::State::cancelled,
                "running queued solve did not report cancellation");
        require(queue.find(cancelledId)->status==dynsyslab::core::JobStatus::cancelled,
                "cancelled status was overwritten");
        const auto workbenchCancelledId=queue.submit(solve(1000.0,0.001));
        const auto workbenchCancelled=queue.runNext("workbench-cancel-test",{}, {},[]{return true;});
        require(workbenchCancelled.state==dynsyslab::workflow::QueueRunResult::State::cancelled
                &&queue.find(workbenchCancelledId)->status==dynsyslab::core::JobStatus::cancelled,
                "external workbench cancellation did not persist queue state");
        bool rejected=false;try{(void)queue.cancel(completedId);}
        catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"terminal job cancellation must be rejected");

        std::filesystem::remove_all(root,ignored);
        std::cout<<"dynsyslab_workflow_queue_tests: PASS\n";return 0;
    } catch(const std::exception& error) {
        std::cerr<<"dynsyslab_workflow_queue_tests: FAIL: "<<error.what()<<'\n';
        std::filesystem::remove_all(root,ignored);return 1;
    }
}
