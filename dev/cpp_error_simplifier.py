import re
from collections import defaultdict


def simplify_cpp_error(error_log: str) -> str:
    simplified_output = defaultdict(lambda: defaultdict(set))

    # Define various regex patterns
    filepath_pattern = re.compile(r'^(.*\.cpp|.*\.hpp|.*\.h|.*\.cc):\d+:\d+:')
    error_type_pattern = re.compile(r'error:|warning:')
    instantiation_pattern = re.compile(r'in instantiation of')

    current_file = None  # Initialize key to None

    for line in error_log.strip().split("\n"):
        # Match file paths and line numbers
        filepath_match = filepath_pattern.match(line)
        if filepath_match:
            current_file = filepath_match.group(0)
            continue

        # Skip if current_file is not yet set
        if current_file is None:
            continue

        # Match error or warning types
        error_type_match = error_type_pattern.search(line)
        if error_type_match:
            simplified_output[current_file]['Type'].add(error_type_match.group(0))
            continue

        # Match instantiation related errors
        instantiation_match = instantiation_pattern.search(line)
        if instantiation_match:
            simplified_output[current_file]['Instantiation'].add(line.strip())
            continue

        # Capture other lines that may contain useful information
        simplified_output[current_file]['Misc'].add(line.strip())

    # Generate the final simplified output
    final_output = []
    for file_key, details in simplified_output.items():
        final_output.append(f"File: {file_key}")
        for detail_type, detail_lines in details.items():
            final_output.append(f"{detail_type}: {', '.join(detail_lines)}")

    return "\n".join(final_output)


if __name__ == "__main__":
    error_log = """
ERROR: /Users/geoffreygarrett/CLionProjects/bor/odin/examples/BUILD:14:10: Compiling odin/examples/example_dev_autodiff.cpp failed: (Exit 1): sandbox-exec failed: error executing command 
  (cd /private/var/tmp/_bazel_geoffreygarrett/c4c7b397206bd8d6eaedc6c5e80968d7/sandbox/darwin-sandbox/122/execroot/bor && \
  exec env - \
    APPLE_SDK_PLATFORM=MacOSX \
    APPLE_SDK_VERSION_OVERRIDE=14.0 \
    ARCHFLAGS='-arch arm64' \
    DEVELOPER_DIR=/Applications/Xcode-beta.app/Contents/Developer \
    PATH=/Users/geoffreygarrett/mambaforge/bin:/Users/geoffreygarrett/mambaforge/condabin:/Users/geoffreygarrett/google-cloud-sdk/bin:/opt/homebrew/bin:/opt/homebrew/sbin:/usr/local/bin:/System/Cryptexes/App/usr/bin:/usr/bin:/bin:/usr/sbin:/sbin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/local/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/appleinternal/bin:/Library/Apple/usr/bin \
    SDKROOT=/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk \
    TMPDIR=/var/folders/mr/jfbmc_8j7r93fs_907ycq9vm0000gn/T/ \
    XCODE_VERSION_OVERRIDE=15.0.0.15A5229h \
    ZERO_AR_DATE=1 \
  /usr/bin/sandbox-exec -f /private/var/tmp/_bazel_geoffreygarrett/c4c7b397206bd8d6eaedc6c5e80968d7/sandbox/darwin-sandbox/122/sandbox.sb /var/tmp/_bazel_geoffreygarrett/install/da3a4a03da022081a13fc9b7c75058bd/process-wrapper '--timeout=0' '--kill_delay=15' '--stats=/private/var/tmp/_bazel_geoffreygarrett/c4c7b397206bd8d6eaedc6c5e80968d7/sandbox/darwin-sandbox/122/stats.out' external/local_config_apple_cc/wrapped_clang_pp '-D_FORTIFY_SOURCE=1' -fstack-protector -fcolor-diagnostics -Wall -Wthread-safety -Wself-assign -fno-omit-frame-pointer -g0 -O2 -DNDEBUG '-DNS_BLOCK_ASSERTIONS=1' '-std=c++14' '-fdebug-prefix-map=__BAZEL_EXECUTION_ROOT__=.' '-fdebug-prefix-map=__BAZEL_XCODE_DEVELOPER_DIR__=/PLACEHOLDER_DEVELOPER_DIR' -iquote . -iquote bazel-out/darwin_arm64-opt/bin -iquote external/com_github_eigen_eigen -iquote bazel-out/darwin_arm64-opt/bin/external/com_github_eigen_eigen -iquote external/com_github_gabime_spdlog -iquote bazel-out/darwin_arm64-opt/bin/external/com_github_gabime_spdlog -iquote external/com_github_fmtlib_fmt -iquote bazel-out/darwin_arm64-opt/bin/external/com_github_fmtlib_fmt -iquote external/local_brew_gnu_gsl -iquote bazel-out/darwin_arm64-opt/bin/external/local_brew_gnu_gsl -iquote external/com_github_oneapi_onetbb -iquote bazel-out/darwin_arm64-opt/bin/external/com_github_oneapi_onetbb -iquote external/com_github_uscilab_cereal -iquote bazel-out/darwin_arm64-opt/bin/external/com_github_uscilab_cereal -iquote external/com_github_autodiff_autodiff -iquote bazel-out/darwin_arm64-opt/bin/external/com_github_autodiff_autodiff -iquote external/bazel_tools -iquote bazel-out/darwin_arm64-opt/bin/external/bazel_tools -isystem odin/include -isystem bazel-out/darwin_arm64-opt/bin/odin/include -isystem external/com_github_eigen_eigen -isystem bazel-out/darwin_arm64-opt/bin/external/com_github_eigen_eigen -isystem external/com_github_gabime_spdlog/include -isystem bazel-out/darwin_arm64-opt/bin/external/com_github_gabime_spdlog/include -isystem external/com_github_fmtlib_fmt/include -isystem bazel-out/darwin_arm64-opt/bin/external/com_github_fmtlib_fmt/include -isystem external/local_brew_gnu_gsl/include -isystem bazel-out/darwin_arm64-opt/bin/external/local_brew_gnu_gsl/include -isystem external/com_github_oneapi_onetbb/include -isystem bazel-out/darwin_arm64-opt/bin/external/com_github_oneapi_onetbb/include -isystem external/com_github_uscilab_cereal/include -isystem bazel-out/darwin_arm64-opt/bin/external/com_github_uscilab_cereal/include -isystem external/com_github_autodiff_autodiff -isystem bazel-out/darwin_arm64-opt/bin/external/com_github_autodiff_autodiff -MD -MF bazel-out/darwin_arm64-opt/bin/odin/examples/_objs/example_dev_autodiff/example_dev_autodiff.d -DFMT_HEADER_ONLY -DSPDLOG_FMT_EXTERNAL -DUSE_PTHREAD -D_XOPEN_SOURCE -DODIN_AUTODIFF '-DBAZEL_CURRENT_REPOSITORY=""' '-frandom-seed=bazel-out/darwin_arm64-opt/bin/odin/examples/_objs/example_dev_autodiff/example_dev_autodiff.o' -isysroot __BAZEL_XCODE_SDKROOT__ -F__BAZEL_XCODE_SDKROOT__/System/Library/Frameworks -F__BAZEL_XCODE_DEVELOPER_DIR__/Platforms/MacOSX.platform/Developer/Library/Frameworks -no-canonical-prefixes -pthread '-std=c++20' -no-canonical-prefixes -Wno-builtin-macro-redefined '-D__DATE__="redacted"' '-D__TIMESTAMP__="redacted"' '-D__TIME__="redacted"' -target arm64-apple-macosx14.0 -c odin/examples/example_dev_autodiff.cpp -o bazel-out/darwin_arm64-opt/bin/odin/examples/_objs/example_dev_autodiff/example_dev_autodiff.o)
In file included from odin/examples/example_dev_autodiff.cpp:3:
In file included from external/com_github_eigen_eigen/Eigen/Core:50:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/complex:243:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/sstream:191:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/istream:165:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/ostream:170:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/bitset:131:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/string:560:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/__memory_resource/polymorphic_allocator.h:20:
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/tuple:1817:5: error: attempt to use a deleted function
    _VSTD::__invoke(
    ^
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/__config:715:17: note: expanded from macro '_VSTD'
#  define _VSTD std
                ^
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/tuple:1814:26: note: in instantiation of exception specification for '__apply_tuple_impl<autodiff::detail::Real<1, double> (*const &)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), const std::tuple<autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &, 0UL, 1UL>' requested here
constexpr decltype(auto) __apply_tuple_impl(_Fn && __f, _Tuple && __t,
                         ^
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/tuple:1826:12: note: in instantiation of function template specialization 'std::__apply_tuple_impl<autodiff::detail::Real<1, double> (*const &)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), const std::tuple<autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &, 0UL, 1UL>' requested here
    _VSTD::__apply_tuple_impl(
           ^
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/tuple:1824:26: note: in instantiation of exception specification for 'apply<autodiff::detail::Real<1, double> (*const &)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), const std::tuple<autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &>' requested here
constexpr decltype(auto) apply(_Fn && __f, _Tuple && __t)
                         ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/derivative.hpp:161:19: note: in instantiation of function template specialization 'std::apply<autodiff::detail::Real<1, double> (*const &)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), const std::tuple<autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &>' requested here
    auto u = std::apply(f, at.args);
                  ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:100:13: note: in instantiation of function template specialization 'autodiff::detail::eval<autodiff::detail::Real<1, double> (*)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &>' requested here
        u = eval(f, at, detail::wrt(xi)); // evaluate u with xi seeded so that du/dxi is also computed
            ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:75:21: note: (skipping 4 contexts in backtrace; use -ftemplate-backtrace-limit=0 to see all)
                    f(i++, item[j]);
                    ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:111:5: note: in instantiation of function template specialization 'autodiff::detail::For<0UL, 3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
    For<0, iend>(std::forward<Function>(f));
    ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:5: note: in instantiation of function template specialization 'autodiff::detail::For<3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
    For<N>([&](auto i) constexpr {
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:5: note: in instantiation of function template specialization 'autodiff::detail::ForEach<const std::tuple<Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &, (lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:23)>' requested here
    ForEach(wrt.args, [&](auto& item) constexpr
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:5: note: in instantiation of function template specialization 'autodiff::detail::ForEachWrtVar<(lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:24), Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &>' requested here
    ForEachWrtVar(wrt, [&](auto&& i, auto&& xi) constexpr
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:114:5: note: in instantiation of function template specialization 'autodiff::detail::gradient<autodiff::detail::Real<1, double> (*)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, Eigen::Matrix<double, -1, 1>>' requested here
    gradient(f, wrt, at, u, g);
    ^
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/__type_traits/nat.h:26:5: note: '~__nat' has been explicitly marked deleted here
    ~__nat() = delete;
    ^
In file included from odin/examples/example_dev_autodiff.cpp:3:
In file included from external/com_github_eigen_eigen/Eigen/Core:50:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/complex:243:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/sstream:191:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/istream:165:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/ostream:170:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/bitset:131:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/string:560:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/__memory_resource/polymorphic_allocator.h:20:
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/tuple:1817:5: error: attempt to use a deleted function
    _VSTD::__invoke(
    ^
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/__config:715:17: note: expanded from macro '_VSTD'
#  define _VSTD std
                ^
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/tuple:1826:12: note: in instantiation of function template specialization 'std::__apply_tuple_impl<autodiff::detail::Real<1, double> (*const &)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), const std::tuple<autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &, 0UL, 1UL>' requested here
    _VSTD::__apply_tuple_impl(
           ^
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/tuple:1824:26: note: in instantiation of exception specification for 'apply<autodiff::detail::Real<1, double> (*const &)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), const std::tuple<autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &>' requested here
constexpr decltype(auto) apply(_Fn && __f, _Tuple && __t)
                         ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/derivative.hpp:161:19: note: in instantiation of function template specialization 'std::apply<autodiff::detail::Real<1, double> (*const &)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), const std::tuple<autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &>' requested here
    auto u = std::apply(f, at.args);
                  ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:100:13: note: in instantiation of function template specialization 'autodiff::detail::eval<autodiff::detail::Real<1, double> (*)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &>' requested here
        u = eval(f, at, detail::wrt(xi)); // evaluate u with xi seeded so that du/dxi is also computed
            ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:75:21: note: in instantiation of function template specialization 'autodiff::detail::gradient(autodiff::detail::Real<1, double> (*const &)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), const Wrt<Matrix<Real<1, double>, 3, 1, 0, 3, 1> &, Real<1, double> &, Matrix<Real<1, double>, 3, 1, 0, 3, 1> &> &, const At<Real<1, double> &, Matrix<Real<1, double>, 3, 1, 0, 3, 1> &> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<double, -1, 1> &)::(anonymous class)::operator()<int, autodiff::detail::Real<1, double> &>' requested here
                    f(i++, item[j]);
                    ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:141:9: note: (skipping 3 contexts in backtrace; use -ftemplate-backtrace-limit=0 to see all)
        f(std::get<i>(tuple));
        ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:111:5: note: in instantiation of function template specialization 'autodiff::detail::For<0UL, 3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
    For<0, iend>(std::forward<Function>(f));
    ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:5: note: in instantiation of function template specialization 'autodiff::detail::For<3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
    For<N>([&](auto i) constexpr {
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:5: note: in instantiation of function template specialization 'autodiff::detail::ForEach<const std::tuple<Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &, (lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:23)>' requested here
    ForEach(wrt.args, [&](auto& item) constexpr
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:5: note: in instantiation of function template specialization 'autodiff::detail::ForEachWrtVar<(lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:24), Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &>' requested here
    ForEachWrtVar(wrt, [&](auto&& i, auto&& xi) constexpr
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:114:5: note: in instantiation of function template specialization 'autodiff::detail::gradient<autodiff::detail::Real<1, double> (*)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, Eigen::Matrix<double, -1, 1>>' requested here
    gradient(f, wrt, at, u, g);
    ^
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/__type_traits/nat.h:26:5: note: '~__nat' has been explicitly marked deleted here
    ~__nat() = delete;
    ^
In file included from odin/examples/example_dev_autodiff.cpp:3:
In file included from external/com_github_eigen_eigen/Eigen/Core:50:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/complex:243:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/sstream:191:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/istream:165:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/ostream:170:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/bitset:131:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/string:560:
In file included from /Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/__memory_resource/polymorphic_allocator.h:20:
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/tuple:1826:5: error: no matching function for call to '__apply_tuple_impl'
    _VSTD::__apply_tuple_impl(
    ^~~~~~~~~~~~~~~~~~~~~~~~~~
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/__config:715:17: note: expanded from macro '_VSTD'
#  define _VSTD std
                ^
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/tuple:1810:79: note: expanded from macro '_LIBCPP_NOEXCEPT_RETURN'
#define _LIBCPP_NOEXCEPT_RETURN(...) noexcept(noexcept(__VA_ARGS__)) { return __VA_ARGS__; }
                                                                              ^~~~~~~~~~~
external/com_github_autodiff_autodiff/autodiff/forward/utils/derivative.hpp:161:19: note: in instantiation of function template specialization 'std::apply<autodiff::detail::Real<1, double> (*const &)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), const std::tuple<autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &>' requested here
    auto u = std::apply(f, at.args);
                  ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:100:13: note: in instantiation of function template specialization 'autodiff::detail::eval<autodiff::detail::Real<1, double> (*)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &>' requested here
        u = eval(f, at, detail::wrt(xi)); // evaluate u with xi seeded so that du/dxi is also computed
            ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:75:21: note: in instantiation of function template specialization 'autodiff::detail::gradient(autodiff::detail::Real<1, double> (*const &)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), const Wrt<Matrix<Real<1, double>, 3, 1, 0, 3, 1> &, Real<1, double> &, Matrix<Real<1, double>, 3, 1, 0, 3, 1> &> &, const At<Real<1, double> &, Matrix<Real<1, double>, 3, 1, 0, 3, 1> &> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<double, -1, 1> &)::(anonymous class)::operator()<int, autodiff::detail::Real<1, double> &>' requested here
                    f(i++, item[j]);
                    ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:141:9: note: in instantiation of function template specialization 'autodiff::detail::ForEachWrtVar(const Wrt<Matrix<Real<1, double>, 3, 1, 0, 3, 1> &, Real<1, double> &, Matrix<Real<1, double>, 3, 1, 0, 3, 1> &> &, (lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:24) &&)::(anonymous class)::operator()<Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1>>' requested here
        f(std::get<i>(tuple));
        ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:97:9: note: in instantiation of function template specialization 'autodiff::detail::ForEach(const std::tuple<Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &, (lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:23) &&)::(anonymous class)::operator()<autodiff::detail::Index<0>>' requested here
        f(Index<i>{});
        ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:105:5: note: (skipping 1 context in backtrace; use -ftemplate-backtrace-limit=0 to see all)
    AuxFor<ibegin, ibegin, iend>(std::forward<Function>(f));
    ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:111:5: note: in instantiation of function template specialization 'autodiff::detail::For<0UL, 3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
    For<0, iend>(std::forward<Function>(f));
    ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:5: note: in instantiation of function template specialization 'autodiff::detail::For<3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
    For<N>([&](auto i) constexpr {
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:5: note: in instantiation of function template specialization 'autodiff::detail::ForEach<const std::tuple<Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &, (lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:23)>' requested here
    ForEach(wrt.args, [&](auto& item) constexpr
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:5: note: in instantiation of function template specialization 'autodiff::detail::ForEachWrtVar<(lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:24), Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &>' requested here
    ForEachWrtVar(wrt, [&](auto&& i, auto&& xi) constexpr
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:114:5: note: in instantiation of function template specialization 'autodiff::detail::gradient<autodiff::detail::Real<1, double> (*)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, Eigen::Matrix<double, -1, 1>>' requested here
    gradient(f, wrt, at, u, g);
    ^
/Applications/Xcode-beta.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX14.0.sdk/usr/include/c++/v1/tuple:1814:26: note: candidate template ignored: substitution failure [with _Fn = autodiff::detail::Real<1, double> (*const &)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), _Tuple = const std::tuple<autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &, _Id = <0, 1>]
constexpr decltype(auto) __apply_tuple_impl(_Fn && __f, _Tuple && __t,
                         ^
In file included from odin/examples/example_dev_autodiff.cpp:10:
In file included from external/com_github_autodiff_autodiff/autodiff/forward/dual/eigen.hpp:37:
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:80:14: error: no matching function for call to object of type '(lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:24)'
        else f(i++, item); // call given f with current index and variable from item (a number, not a vector)
             ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:141:9: note: in instantiation of function template specialization 'autodiff::detail::ForEachWrtVar(const Wrt<Matrix<Real<1, double>, 3, 1, 0, 3, 1> &, Real<1, double> &, Matrix<Real<1, double>, 3, 1, 0, 3, 1> &> &, (lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:24) &&)::(anonymous class)::operator()<autodiff::detail::Real<1, double>>' requested here
        f(std::get<i>(tuple));
        ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:97:9: note: in instantiation of function template specialization 'autodiff::detail::ForEach(const std::tuple<Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &, (lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:23) &&)::(anonymous class)::operator()<autodiff::detail::Index<1>>' requested here
        f(Index<i>{});
        ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:98:9: note: in instantiation of function template specialization 'autodiff::detail::AuxFor<1UL, 0UL, 3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
        AuxFor<i + 1, ibegin, iend>(std::forward<Function>(f));
        ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:105:5: note: in instantiation of function template specialization 'autodiff::detail::AuxFor<0UL, 0UL, 3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
    AuxFor<ibegin, ibegin, iend>(std::forward<Function>(f));
    ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:111:5: note: in instantiation of function template specialization 'autodiff::detail::For<0UL, 3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
    For<0, iend>(std::forward<Function>(f));
    ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:5: note: in instantiation of function template specialization 'autodiff::detail::For<3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
    For<N>([&](auto i) constexpr {
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:5: note: in instantiation of function template specialization 'autodiff::detail::ForEach<const std::tuple<Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &, (lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:23)>' requested here
    ForEach(wrt.args, [&](auto& item) constexpr
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:5: note: in instantiation of function template specialization 'autodiff::detail::ForEachWrtVar<(lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:24), Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &>' requested here
    ForEachWrtVar(wrt, [&](auto&& i, auto&& xi) constexpr
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:114:5: note: in instantiation of function template specialization 'autodiff::detail::gradient<autodiff::detail::Real<1, double> (*)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, Eigen::Matrix<double, -1, 1>>' requested here
    gradient(f, wrt, at, u, g);
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:24: note: candidate template ignored: substitution failure [with i:auto = int, xi:auto = autodiff::detail::Real<1, double> &]
    ForEachWrtVar(wrt, [&](auto&& i, auto&& xi) constexpr
                       ^
In file included from odin/examples/example_dev_autodiff.cpp:9:
In file included from external/com_github_autodiff_autodiff/autodiff/forward/dual.hpp:33:
In file included from external/com_github_autodiff_autodiff/autodiff/forward/dual/dual.hpp:42:
In file included from external/com_github_autodiff_autodiff/autodiff/common/numbertraits.hpp:33:
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:141:9: error: no matching function for call to object of type '(lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:23)'
        f(std::get<i>(tuple));
        ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:97:9: note: in instantiation of function template specialization 'autodiff::detail::ForEach(const std::tuple<Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &, (lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:23) &&)::(anonymous class)::operator()<autodiff::detail::Index<2>>' requested here
        f(Index<i>{});
        ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:98:9: note: in instantiation of function template specialization 'autodiff::detail::AuxFor<2UL, 0UL, 3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
        AuxFor<i + 1, ibegin, iend>(std::forward<Function>(f));
        ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:98:9: note: in instantiation of function template specialization 'autodiff::detail::AuxFor<1UL, 0UL, 3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:105:5: note: in instantiation of function template specialization 'autodiff::detail::AuxFor<0UL, 0UL, 3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
    AuxFor<ibegin, ibegin, iend>(std::forward<Function>(f));
    ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:111:5: note: in instantiation of function template specialization 'autodiff::detail::For<0UL, 3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
    For<0, iend>(std::forward<Function>(f));
    ^
external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:5: note: in instantiation of function template specialization 'autodiff::detail::For<3UL, (lambda at external/com_github_autodiff_autodiff/autodiff/common/meta.hpp:140:12)>' requested here
    For<N>([&](auto i) constexpr {
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:5: note: in instantiation of function template specialization 'autodiff::detail::ForEach<const std::tuple<Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &> &, (lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:23)>' requested here
    ForEach(wrt.args, [&](auto& item) constexpr
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:5: note: in instantiation of function template specialization 'autodiff::detail::ForEachWrtVar<(lambda at external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:97:24), Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &>' requested here
    ForEachWrtVar(wrt, [&](auto&& i, auto&& xi) constexpr
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:114:5: note: in instantiation of function template specialization 'autodiff::detail::gradient<autodiff::detail::Real<1, double> (*)(Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, const Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &), Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double> &, Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1> &, autodiff::detail::Real<1, double>, Eigen::Matrix<double, -1, 1>>' requested here
    gradient(f, wrt, at, u, g);
    ^
external/com_github_autodiff_autodiff/autodiff/forward/utils/gradient.hpp:67:23: note: candidate template ignored: substitution failure [with item:auto = Eigen::Matrix<autodiff::detail::Real<1, double>, 3, 1, 0, 3, 1>]
    ForEach(wrt.args, [&](auto& item) constexpr
                      ^
5 errors generated.
Error in child process '/usr/bin/xcrun'. 1
Target //odin/examples:example_dev_autodiff failed to build
Use --verbose_failures to see the command lines of failed build steps.
INFO: Elapsed time: 0.917s, Critical Path: 0.77s
INFO: 2 processes: 2 internal.
FAILED: Build did NOT complete successfully
ERROR: Build failed. Not running target
    """
    print(simplify_cpp_error(error_log))
