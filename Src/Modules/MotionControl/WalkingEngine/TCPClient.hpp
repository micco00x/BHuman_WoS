#pragma once

// STL
#include <functional>
#include <thread>

// Unix socket
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <fcntl.h>
#include <unistd.h>
#include <poll.h>

// Eigen
#include <Eigen/Core>

#include "Foot.hpp"
#include "types.hpp"

class TCPClient {
 public:
  bool connectToServer(const char* ip_addr, in_port_t port) {
    sockfd_ = socket(AF_INET, SOCK_STREAM, 0);

    struct sockaddr_in client;
    client.sin_family = AF_INET;
    client.sin_port = port;
    client.sin_addr.s_addr = inet_addr(ip_addr);

    if (connect(sockfd_, (sockaddr*) &client, sizeof(sockaddr_in)) == -1) {
      return false;
    }

    return true;
  }

  template <class F, class... Args>
  void subscribeToFootstepPlan(F&& f, Args&&... args) {
    footstep_plan_thread_ = std::thread(&TCPClient::footstepPlanPoll, this);
    footstepPlanCallback_ = std::bind(f, args..., std::placeholders::_1);
  }

 private:
  
  void footstepPlanPoll() {
    struct pollfd pfds[1];
    pfds[0].fd = sockfd_;
    pfds[0].events = POLLIN;

    while (1) {
      // Poll and receive footstep plan if POLLIN has been set:
      if (poll(pfds, 1, -1) > 0) {
        if (pfds[0].revents & POLLIN) {
          FootstepPlan footstepPlan;
          if (recvFootstepPlan(footstepPlan)) {
            footstepPlanCallback_(footstepPlan);
          } else {
            break;
          }
        }
      } else {
        break;
      }
    }
  }

  bool recvFootstepPlan(FootstepPlan& footstepPlan) {
    ssize_t ss;
    uint32_t footstep_plan_size;

    // Receive size of the footstep plan:
    ss = recv(sockfd_, &footstep_plan_size, sizeof(footstep_plan_size),
        MSG_DONTWAIT);
    if (ss <= 0) {
      return false;
    }

    const size_t double_support_configuration_size =
        sizeof(double) * 9 + sizeof(Foot);
    const size_t buffer_plan_size =
        footstep_plan_size * double_support_configuration_size;
    std::unique_ptr<char[]> buffer_plan(new char[buffer_plan_size]);

    // Receive footstep plan:
    ss = recv(sockfd_, &buffer_plan[0], buffer_plan_size, MSG_WAITALL);
    if (ss <= 0) {
      return false;
    }

    for (size_t k = 0; k < footstep_plan_size; ++k) {
      double q[8], h_z;
      Foot support_foot;
      char* ptr = &buffer_plan[0] + k * double_support_configuration_size;
      memcpy(q, ptr, sizeof(q));
      ptr += sizeof(q);
      memcpy(&support_foot, ptr, sizeof(support_foot));
      ptr += sizeof(support_foot);
      memcpy(&h_z, ptr, sizeof(h_z));
      Eigen::Vector4d qL;
      Eigen::Vector4d qR;
      qL << q[0], q[1], q[2], q[3];
      qR << q[4], q[5], q[6], q[7];
      footstepPlan.push_back(Configuration(qL, qR, support_foot, h_z));
    }

    return true;
  }

  int sockfd_;

  std::thread footstep_plan_thread_;
  std::function<void(const FootstepPlan&)> footstepPlanCallback_;

}; // end class TCPClient
