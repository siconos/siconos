#pragma once

namespace siconos::storage::database {

  template<typename Info, typename UnderlyingDb>
  struct database {
    using info = Info;

    UnderlyingDb _data;
  };

}
